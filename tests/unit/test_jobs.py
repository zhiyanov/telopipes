"""Scheduler state transitions, with execution stubbed out."""
import pytest

from telomers import db, jobs


@pytest.fixture
def scheduler(tmp_path, monkeypatch):
    from telomers.config import settings
    monkeypatch.setattr(settings, "data_root", tmp_path)
    monkeypatch.setattr(db, "_engine", None)
    with db.session() as session:
        session.add(db.Run(id="r1", pipeline="telonp", status="created"))
        session.commit()
    return jobs.Scheduler()


def _status(run_id: str) -> str:
    with db.session() as session:
        return session.get(db.Run, run_id).status


def test_recover_marks_orphaned_runs_interrupted(scheduler):
    """A run still marked running when the process died cannot be adopted. The work is not
    lost: Snakemake's state is in the run directory, so resuming re-executes only the rest."""
    with db.session() as session:
        run = session.get(db.Run, "r1")
        run.status = "running"
        session.add(run)
        session.commit()
    scheduler.recover()
    assert _status("r1") == "interrupted"


def test_cancelling_is_sticky(scheduler):
    """Stopping the container makes it exit non-zero. That must not be reported back to the
    user as a failure when they asked for the cancellation."""
    jobs.Scheduler._set_status("r1", "cancelled")
    jobs.Scheduler._set_status("r1", "running")
    assert _status("r1") == "cancelled"
    jobs.Scheduler._finish("r1", 137, None)
    assert _status("r1") == "cancelled"


def test_resubmission_clears_the_cancelled_state(scheduler, tmp_path):
    """The stickiness above must not make a cancelled run impossible to start again."""
    jobs.Scheduler._set_status("r1", "cancelled")
    scheduler.submit(jobs.Job("r1", "telonp", tmp_path / "reads.fastq", None))
    assert _status("r1") == "queued"


def test_finish_records_success_and_summary(scheduler):
    jobs.Scheduler._finish("r1", 0, {"pipeline": "telonp"})
    with db.session() as session:
        run = session.get(db.Run, "r1")
    assert run.status == "succeeded" and run.exit_code == 0 and run.summary["pipeline"] == "telonp"


def test_finish_records_failure(scheduler):
    jobs.Scheduler._finish("r1", 1, None)
    assert _status("r1") == "failed"


def test_queued_run_is_cancelled_without_touching_podman(scheduler, tmp_path):
    import asyncio
    scheduler.submit(jobs.Job("r1", "telonp", tmp_path / "reads.fastq", None))
    assert asyncio.run(scheduler.cancel("r1")) is True
    assert _status("r1") == "cancelled"
