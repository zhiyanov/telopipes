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

from . import doctor, execute, pipelines
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
    print(f"outputs: {run.work}")
    print(f"log:     {run.log}")
    if code != 0:
        print("re-run the same command to resume: Snakemake skips completed steps.")
    return code


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

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
