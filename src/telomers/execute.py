"""Job execution: write a config, start one container, let Snakemake do the rest.

Snakemake runs *inside* the pipeline image rather than on the host. Every method in and out of
scope is a single-image pipeline -- telomere-r bundles minimap2, samtools and R together -- so
nothing needs an orchestrator that spans images. That keeps Snakemake and its ~40 dependencies
off the host, and means this module makes exactly one `podman run` call per run: no socket, no
nested containers, no host/container path mismatch. See docs/PLAN.md sections 2 and 3.
"""
from __future__ import annotations

import subprocess
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import yaml
from pydantic import BaseModel

from .config import settings
from .pipelines import Pipeline

#: Fixed paths inside the container. Inputs are mounted read-only as individual files, and the
#: run's own work/ is the only writable mount -- never the repo, never the data root.
READS_IN = "/inputs/reads.fastq"
REFERENCE_IN = "/refs/reference.fasta"
WORK_IN = "/work"


@dataclass(frozen=True)
class Run:
    id: str
    pipeline: Pipeline
    directory: Path

    @property
    def work(self) -> Path:
        return self.directory / "work"

    @property
    def log(self) -> Path:
        return self.directory / "run.log"

    @property
    def config(self) -> Path:
        return self.work / "config.yaml"


def new_run_id() -> str:
    return f"{datetime.now(timezone.utc):%Y%m%d-%H%M%S}-{uuid.uuid4().hex[:6]}"


def prepare(
    pipeline: Pipeline,
    params: BaseModel,
    reads: Path,
    reference: Path | None = None,
    run_id: str | None = None,
) -> Run:
    """Create the run directory and write the config Snakemake will read."""
    if pipeline.requires_reference and reference is None:
        raise ValueError(f"pipeline {pipeline.id} requires a reference genome")

    rid = run_id or new_run_id()
    run = Run(rid, pipeline, settings.runs_dir() / rid)
    run.work.mkdir(parents=True, exist_ok=True)

    config = params.model_dump()
    config["reads"] = READS_IN
    if reference is not None:
        config["reference"] = REFERENCE_IN
    run.config.write_text(yaml.safe_dump(config, sort_keys=True))

    # Record what produced this run, next to it, so the run directory is self-describing and
    # the database can be rebuilt by walking the data root.
    (run.directory / "run.json").write_text(
        yaml.safe_dump(
            {
                "run_id": run.id,
                "pipeline": pipeline.id,
                "pipeline_version": pipeline.version,
                "image": pipeline.image,
                "created": datetime.now(timezone.utc).isoformat(),
                "inputs": {
                    "reads": str(reads.resolve()),
                    "reference": str(reference.resolve()) if reference else None,
                },
                "params": config,
            },
            sort_keys=True,
        )
    )
    return run


def build_command(
    run: Run, reads: Path, reference: Path | None, *, cores: int, memory: str
) -> list[str]:
    """The exact podman invocation. Kept separate so it can be asserted on in tests."""
    cmd = [settings.container_engine, "run", "--rm", "--name", f"telomers-{run.id}"]
    if settings.platform:
        cmd += ["--platform", settings.platform]
    cmd += ["--memory", memory, "--cpus", str(cores)]
    # ,z is required: SELinux is enforcing here, and without it the mount surfaces as
    # "Permission denied" inside the container, which reads like a pipeline bug.
    cmd += ["-v", f"{run.work}:{WORK_IN}:z"]
    cmd += ["-v", f"{reads.resolve()}:{READS_IN}:ro,z"]
    if reference is not None:
        cmd += ["-v", f"{reference.resolve()}:{REFERENCE_IN}:ro,z"]
    cmd += ["-w", WORK_IN, run.pipeline.image]
    cmd += [
        "snakemake",
        "--snakefile", "/opt/telomers/Snakefile",
        "--cores", str(cores),
        "--rerun-incomplete",
        "--nolock",
    ]
    return cmd


def execute(
    run: Run,
    reads: Path,
    reference: Path | None,
    *,
    cores: int | None = None,
    memory: str | None = None,
    echo: bool = True,
) -> int:
    """Run to completion, teeing output to run.log. Returns the exit code.

    Re-invoking this on an existing run directory resumes it: Snakemake skips outputs that
    already exist, so an interrupted pipeline-A run picks up at `analyze` without redoing the
    alignment. That is the whole of our resume story.
    """
    cmd = build_command(
        run, reads, reference,
        cores=cores or settings.threads,
        memory=memory or settings.memory,
    )
    if echo:
        print("$ " + " ".join(cmd), flush=True)

    with run.log.open("a") as logfile:
        logfile.write(f"\n$ {' '.join(cmd)}\n")
        logfile.flush()
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
        )
        assert process.stdout is not None
        for line in process.stdout:
            logfile.write(line)
            logfile.flush()
            if echo:
                print(line, end="", flush=True)
        return process.wait()


def cancel(run_id: str) -> None:
    """Stop a running container. One container per run, so one name to stop."""
    subprocess.run(
        [settings.container_engine, "stop", "--time", "10", f"telomers-{run_id}"],
        check=False, capture_output=True,
    )


# --------------------------------------------------------------------------- reference prep

def cut_reference(
    genome: Path,
    out_fasta: Path,
    *,
    arm_length: int = 500_000,
    haplotype: str = "all",
    pipeline_id: str = "telomere-r",
    echo: bool = True,
) -> str:
    """Derive a cut reference from a whole genome, in the pipeline's own container.

    Run in the image rather than on the host so the cutting is done by exactly the code that
    ships with the pipeline -- the same file, not a copy that can drift. It is pure standard
    library, so the container costs nothing but the round trip.
    """
    from . import pipelines as _pipelines

    pipeline = _pipelines.get(pipeline_id)
    genome = genome.resolve()
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_dir = out_fasta.parent.resolve()

    cmd = [settings.container_engine, "run", "--rm"]
    if settings.platform:
        cmd += ["--platform", settings.platform]
    cmd += [
        "-v", f"{genome}:/genome/{genome.name}:ro,z",
        "-v", f"{out_dir}:/out:z",
        "-w", "/out", pipeline.image,
        "python3", "/opt/telomers/cut_reference.py",
        f"/genome/{genome.name}", f"/out/{out_fasta.name}",
        "--arm-length", str(arm_length),
        "--haplotype", haplotype,
        "--report", f"/out/{out_fasta.stem}.report.json",
    ]
    if echo:
        print("$ " + " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if echo and proc.stdout:
        print(proc.stdout, end="")
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or proc.stdout or "cutting failed").strip())
    return proc.stdout
