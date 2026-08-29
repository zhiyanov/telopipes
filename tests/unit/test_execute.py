from pathlib import Path

from telomers import execute, pipelines
from telomers.params import TeloNPParams


def _run(tmp_path):
    pipeline = pipelines.get("telonp")
    run = execute.Run("test-run", pipeline, tmp_path / "test-run")
    run.work.mkdir(parents=True)
    return run


def test_command_mounts_inputs_read_only_and_labels_them(tmp_path):
    run = _run(tmp_path)
    reads = tmp_path / "reads.fastq"
    reads.touch()
    cmd = execute.build_command(run, reads, None, cores=4, memory="8g")
    joined = " ".join(cmd)

    # Inputs read-only, the run's own work/ writable, nothing else mounted.
    assert f"{reads.resolve()}:{execute.READS_IN}:ro,z" in joined
    assert f"{run.work}:{execute.WORK_IN}:z" in joined
    assert sum(1 for a in cmd if a == "-v") == 2

    # ,z on every mount: SELinux is enforcing, and a missing label shows up as a
    # permission error inside the container that reads like a pipeline bug.
    mounts = [cmd[i + 1] for i, a in enumerate(cmd) if a == "-v"]
    assert mounts and all(m.endswith(",z") or m.endswith(":z") for m in mounts)

    # No host path leaks into the command Snakemake runs.
    assert cmd[-4:] == ["--cores", "4", "--rerun-incomplete", "--nolock"]
    assert str(tmp_path) not in " ".join(cmd[cmd.index(run.pipeline.image):])


def test_reference_is_mounted_only_when_given(tmp_path):
    run = _run(tmp_path)
    reads = tmp_path / "reads.fastq"
    ref = tmp_path / "ref.fasta"
    reads.touch()
    ref.touch()
    assert execute.REFERENCE_IN not in " ".join(execute.build_command(run, reads, None, cores=1, memory="1g"))
    assert execute.REFERENCE_IN in " ".join(execute.build_command(run, reads, ref, cores=1, memory="1g"))


def test_prepare_writes_container_paths_not_host_paths(tmp_path, monkeypatch):
    monkeypatch.setattr(execute.settings, "data_root", tmp_path)
    reads = tmp_path / "reads.fastq"
    reads.touch()
    run = execute.prepare(pipelines.get("telonp"), TeloNPParams(sample_label="s"), reads,
                          run_id="r1")
    config = run.config.read_text()
    # Snakemake resolves these inside the container, so they must be the container's paths.
    assert execute.READS_IN in config
    assert str(tmp_path) not in config


def test_prepare_rejects_a_missing_reference(tmp_path, monkeypatch):
    monkeypatch.setattr(execute.settings, "data_root", tmp_path)
    reads = tmp_path / "reads.fastq"
    reads.touch()
    try:
        execute.prepare(pipelines.get("telomere-r"),
                        TeloNPParams(sample_label="s"), reads, None)
    except ValueError as exc:
        assert "requires a reference" in str(exc)
    else:
        raise AssertionError("expected ValueError")
