"""Command line interface.

`telomers run` replaces the reference project's run_all.sh: instead of editing a shell script
to change sample, barcode or thresholds, they are arguments.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from pydantic import ValidationError

from . import db, doctor, execute, pipelines, plots, results, storage, validate
from .config import settings


def _add_params(parser: argparse.ArgumentParser, models) -> None:
    """Expose the Pydantic params models as CLI flags, so the models stay the only definition.

    Flags are the union across pipelines (several are shared, e.g. sample_label), and none is
    marked required here -- validation belongs to the model, which knows which pipeline is
    actually being run and can say so properly.
    """
    seen: set[str] = set()
    for model in models:
        for name, field in model.model_fields.items():
            if name in seen:
                continue
            seen.add(name)
            default = "required" if field.is_required() else field.default
            parser.add_argument(
                "--" + name.replace("_", "-"),
                default=None,
                help=f"{field.description or name} [{default}]",
            )


def cmd_pipelines(args: argparse.Namespace) -> int:
    for pipeline in sorted(pipelines.registry().values(), key=lambda p: p.id):
        print(f"{pipeline.id}  (v{pipeline.version})")
        print(f"  {pipeline.title}")
        if pipeline.summary:
            for line in pipeline.summary.split("\n"):
                print(f"    {line}")
        print(f"  reference required: {'yes' if pipeline.requires_reference else 'no'}")
        print(f"  produces:           {', '.join(pipeline.provides)}")
        print(f"  image:              {pipeline.image}")
        if pipeline.citation:
            print(f"  citation:           {pipeline.citation}")
        print()
    return 0


def cmd_doctor(args: argparse.Namespace) -> int:
    symbols = {doctor.OK: "  ok  ", doctor.WARN: " warn ", doctor.FAIL: " FAIL "}
    worst = 0
    for check in doctor.run_checks():
        print(f"[{symbols[check.status]}] {check.name:22} {check.detail}")
        if check.status == doctor.FAIL:
            worst = 1
    return worst


def cmd_images(args: argparse.Namespace) -> int:
    reg = pipelines.registry()
    targets = [reg[args.pipeline]] if args.pipeline else sorted(reg.values(), key=lambda p: p.id)
    for pipeline in targets:
        cmd = [settings.container_engine, "build"]
        if settings.platform:
            cmd += ["--platform", settings.platform]
        cmd += ["-t", pipeline.image, str(pipeline.directory)]
        print("$ " + " ".join(cmd), flush=True)
        code = subprocess.call(cmd)
        if code != 0:
            print(f"build failed for {pipeline.id}", file=sys.stderr)
            return code
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    pipeline = pipelines.get(args.pipeline)
    model = pipeline.params_model

    # args carries the union of every pipeline's flags; keep only this model's fields.
    supplied = {
        name: getattr(args, name)
        for name in model.model_fields
        if getattr(args, name, None) is not None
    }
    supplied.setdefault("sample_label", Path(args.reads).name.split(".")[0])

    try:
        params = model(**supplied)
    except ValidationError as exc:
        for err in exc.errors():
            field = ".".join(str(p) for p in err["loc"]) or "(model)"
            print(f"error: {field}: {err['msg']}", file=sys.stderr)
        return 2

    reads = Path(args.reads)
    if not reads.is_file():
        print(f"error: reads not found: {reads}", file=sys.stderr)
        return 2
    reference = Path(args.reference) if args.reference else None
    if reference is not None and not reference.is_file():
        print(f"error: reference not found: {reference}", file=sys.stderr)
        return 2

    try:
        run = execute.prepare(pipeline, params, reads, reference, run_id=args.run_id)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"run {run.id}  ->  {run.directory}")
    code = execute.execute(run, reads, reference,
                           cores=args.cores, memory=args.memory)
    print(f"\nrun {run.id} {'succeeded' if code == 0 else f'FAILED (exit {code})'}")

    if code == 0:
        out = run.directory / "results"
        try:
            normalised = results.normalise(
                run.work, out, pipeline_id=pipeline.id,
                extra={"sample_label": params.sample_label, "run_id": run.id},
            )
            made = plots.render_all(normalised.per_read, out, sample=params.sample_label)
            print(f"canonical results: {out}  ({len(made)} figures)")
            _print_summary(json.loads(normalised.summary.read_text()))
        except Exception as exc:  # noqa: BLE001 - the run itself succeeded; say so and move on
            print(f"warning: could not normalise outputs: {exc}", file=sys.stderr)

    print(f"outputs: {run.work}")
    print(f"log:     {run.log}")
    if code != 0:
        print("re-run the same command to resume: Snakemake skips completed steps.")
    return code


def _print_summary(summary: dict) -> None:
    reads = summary["reads"]
    tel = summary["telomere_bp"]
    print(f"\n  reads      {reads['passed']:,} of {reads['total']:,} passed "
          f"({reads['pass_rate']:.1%})")
    if tel["n"]:
        print(f"  telomere   median {tel['median']:,.1f} bp   mean {tel['mean']:,.1f} bp   "
              f"sd {tel['sd']:,.1f}")
        print(f"             range {tel['min']:,}-{tel['max']:,} bp   "
              f"IQR {tel['q25']:,.0f}-{tel['q75']:,.0f}")
    if summary["qc_breakdown"]:
        drops = "  ".join(f"{k.replace('_', ' ')} {v:,}" for k, v in summary["qc_breakdown"].items())
        print(f"  dropped    {drops}")


def cmd_datasets(args: argparse.Namespace) -> int:
    if args.action == "add":
        source = Path(args.path)
        if args.registered:
            storage.check_import_allowed(source)
        meta = storage.register_dataset(source, name=args.name, link=not args.copy)
        with db.session() as s:
            s.add(db.from_meta(db.Dataset, meta))
            s.commit()
        print(f"{meta['id']}  {meta['name']}")
        print(f"  {meta['n_reads']:,} reads, {meta['total_bases']:,} bases"
              if meta["n_reads"] else "  (not counted)")
        if meta["guessed_barcode"]:
            print(f"  barcode guessed from filename: {meta['guessed_barcode']}")
        print(f"  {'hardlinked' if meta['linked'] else 'copied'} to {meta['path']}")
        return 0

    with db.session() as s:
        rows = s.exec(db.select(db.Dataset)).all()
    if not rows:
        print("no datasets registered")
        return 0
    for row in rows:
        reads = f"{row.n_reads:,} reads" if row.n_reads else "uncounted"
        print(f"{row.id:34} {row.name:36} {reads}")
    return 0


def cmd_references(args: argparse.Namespace) -> int:
    if args.action == "add":
        source = Path(args.path)
        if args.registered:
            storage.check_import_allowed(source)
        try:
            meta = storage.register_reference(source, name=args.name, link=not args.copy)
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        with db.session() as s:
            s.add(db.from_meta(db.Reference, meta))
            s.commit()
        print(f"{meta['id']}  {meta['name']}")
        print(f"  {meta['n_records']} records, arm length {meta['arm_length_bp']:,} bp"
              f"  haplotypes: {meta['haplotypes'] or 'none'}")
        for warning in meta["warnings"]:
            print(f"  warning: {warning}")
        return 0

    if args.action == "check":
        info = validate.validate_reference(Path(args.path))
        print(f"records:    {info.n_records}")
        print(f"arm length: {info.arm_length_bp:,} bp" if info.arm_length_bp
              else "arm length: not uniform")
        print(f"haplotypes: {', '.join(info.haplotypes) or 'none'}")
        for e in info.errors:
            print(f"ERROR: {e}")
        for w in info.warnings:
            print(f"warn:  {w}")
        print("\nusable" if info.ok else "\nNOT usable")
        return 0 if info.ok else 1

    with db.session() as s:
        rows = s.exec(db.select(db.Reference)).all()
    if not rows:
        print("no references registered")
        return 0
    for row in rows:
        print(f"{row.id:34} {row.name:36} {row.n_records:>4} records  "
              f"arm {row.arm_length_bp:,} bp  {row.haplotypes}")
    return 0


def cmd_results(args: argparse.Namespace) -> int:
    directory = storage.runs_dir() / args.run_id
    summary_path = directory / "results" / "summary.json"
    if not summary_path.is_file():
        print(f"no canonical results for run {args.run_id!r} "
              f"(looked in {summary_path})", file=sys.stderr)
        return 2
    summary = json.loads(summary_path.read_text())
    print(f"run {args.run_id}   pipeline {summary['pipeline']}"
          f"   sample {summary.get('sample_label', '?')}")
    _print_summary(summary)

    per_arm = directory / "results" / "per_arm.csv"
    if per_arm.is_file() and args.by_arm:
        import csv as _csv
        with per_arm.open(newline="") as handle:
            rows = list(_csv.DictReader(handle))
        print(f"\n  {'arm':<10}{'n':>6}{'median':>10}{'mean':>10}{'sd':>10}")
        for row in rows:
            print(f"  {row['chrom_arm']:<10}{row['n']:>6}{row['median_telomere_bp']:>10}"
                  f"{row['mean_telomere_bp']:>10}{row['sd_telomere_bp']:>10}")
    elif per_arm.is_file():
        print(f"\n  per-arm results available: --by-arm to show, or {per_arm}")
    return 0


def cmd_runs(args: argparse.Namespace) -> int:
    with db.session() as s:
        rows = s.exec(db.select(db.Run)).all()
    if not rows:
        print("no runs recorded (try: telomers reindex)")
        return 0
    for row in rows:
        summary = row.summary or {}
        tel = summary.get("telomere_bp") or {}
        median = f"median {tel['median']:,.0f} bp" if tel.get("median") else ""
        print(f"{row.id:34} {row.pipeline:12} {row.status:12} {median}")
    return 0


def cmd_reindex(args: argparse.Namespace) -> int:
    counts = db.reindex()
    print("rebuilt index from disk: "
          + ", ".join(f"{v} {k}" for k, v in counts.items()))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="telomers", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("pipelines", help="list available pipelines").set_defaults(func=cmd_pipelines)
    sub.add_parser("doctor", help="check the environment").set_defaults(func=cmd_doctor)

    p_images = sub.add_parser("images", help="build pipeline images")
    p_images.add_argument("action", choices=["build"])
    p_images.add_argument("--pipeline", help="build only this one")
    p_images.set_defaults(func=cmd_images)

    p_run = sub.add_parser("run", help="run a pipeline")
    p_run.add_argument("pipeline", choices=sorted(pipelines.registry()))
    p_run.add_argument("--reads", required=True, help="input FASTQ (optionally gzipped)")
    p_run.add_argument("--reference", help="cut reference FASTA (pipeline A only)")
    p_run.add_argument("--run-id", help="reuse an existing run directory, i.e. resume it")
    p_run.add_argument("--cores", type=int, default=None)
    p_run.add_argument("--memory", default=None)
    # Parameter flags come from the Pydantic models, so they cannot drift from validation.
    _add_params(p_run, {p.params_model for p in pipelines.registry().values()})
    p_run.set_defaults(func=cmd_run)

    p_ds = sub.add_parser("datasets", help="register and list read sets")
    p_ds.add_argument("action", choices=["add", "list"])
    p_ds.add_argument("path", nargs="?", help="FASTQ to register")
    p_ds.add_argument("--name", help="label (default: filename)")
    p_ds.add_argument("--copy", action="store_true",
                      help="copy instead of hardlinking (needed across filesystems)")
    p_ds.add_argument("--registered", action="store_true",
                      help="the path is server-side; enforce the import allowlist")
    p_ds.set_defaults(func=cmd_datasets)

    p_ref = sub.add_parser("references", help="register, check and list cut references")
    p_ref.add_argument("action", choices=["add", "list", "check"])
    p_ref.add_argument("path", nargs="?", help="FASTA to register or check")
    p_ref.add_argument("--name")
    p_ref.add_argument("--copy", action="store_true")
    p_ref.add_argument("--registered", action="store_true")
    p_ref.set_defaults(func=cmd_references)

    p_res = sub.add_parser("results", help="show a run's canonical summary")
    p_res.add_argument("run_id")
    p_res.add_argument("--by-arm", action="store_true", help="also list per-arm statistics")
    p_res.set_defaults(func=cmd_results)

    sub.add_parser("runs", help="list recorded runs").set_defaults(func=cmd_runs)
    sub.add_parser("reindex", help="rebuild the index from disk").set_defaults(func=cmd_reindex)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
