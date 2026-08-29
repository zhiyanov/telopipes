"""Run queue.

One run at a time. Pipeline A's analysis step budgets 10 of this machine's 15 GB, and nobody in
a lab of a few people needs two concurrent runs badly enough to justify the memory pressure or
the code to manage it.

There is no Celery and no Redis: the containers do all the real work, so the event loop only
has to supervise. State lives in the database rather than in memory, which is what makes
recovery after a restart possible -- and the database is itself rebuildable from disk, so a
run's history survives even losing that.
"""
from __future__ import annotations

import asyncio
import json
import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from sqlmodel import select

from . import db, execute, pipelines, plots, results

log = logging.getLogger("telomers.jobs")


@dataclass
class Job:
    run_id: str
    pipeline_id: str
    reads: Path
    reference: Path | None
    cores: int | None = None
    memory: str | None = None


class Scheduler:
    def __init__(self) -> None:
        self._queue: asyncio.Queue[Job] = asyncio.Queue()
        self._worker: asyncio.Task | None = None
        self._current: str | None = None

    # ---------------------------------------------------------------- lifecycle

    async def start(self) -> None:
        self.recover()
        self._worker = asyncio.create_task(self._run_forever())

    async def stop(self) -> None:
        if self._worker:
            self._worker.cancel()
            try:
                await self._worker
            except asyncio.CancelledError:
                pass

    def recover(self) -> None:
        """Anything left running when we died cannot be adopted, so mark it interrupted.

        The work is not lost: Snakemake's state lives in the run's own work directory, so
        resuming re-executes only what did not finish.
        """
        with db.session() as session:
            for run in session.exec(
                select(db.Run).where(db.Run.status.in_(("running", "queued")))
            ).all():
                run.status = "interrupted"
                session.add(run)
            session.commit()

    # ---------------------------------------------------------------- queueing

    def submit(self, job: Job) -> None:
        # force=True because a resubmitted run may currently be "cancelled", and the guard in
        # _set_status exists to stop a dying run from overwriting that -- not to stop the user
        # from starting it again.
        self._set_status(job.run_id, "queued", force=True)
        self._queue.put_nowait(job)

    @property
    def current(self) -> str | None:
        return self._current

    def queued_ahead_of(self, run_id: str) -> int:
        return sum(1 for job in self._queue._queue if job.run_id == run_id)  # noqa: SLF001

    async def cancel(self, run_id: str) -> bool:
        """Cancel a queued or running run.

        A queued run never reaches podman, so it is cancelled in the database alone; a running
        one has exactly one container, named after it.
        """
        if self._current == run_id:
            # Mark it cancelled *before* stopping the container. Stopping makes the process
            # exit non-zero, and _finish would otherwise report that as a failure -- which is
            # what a user cancelling a run least wants to see.
            self._set_status(run_id, "cancelled")
            await asyncio.to_thread(execute.cancel, run_id)
            return True
        remaining = [j for j in self._queue._queue if j.run_id != run_id]  # noqa: SLF001
        if len(remaining) != len(self._queue._queue):  # noqa: SLF001
            self._queue._queue.clear()  # noqa: SLF001
            for job in remaining:
                self._queue._queue.append(job)  # noqa: SLF001
            self._set_status(run_id, "cancelled")
            return True
        return False

    # ---------------------------------------------------------------- worker

    async def _run_forever(self) -> None:
        while True:
            job = await self._queue.get()
            self._current = job.run_id
            try:
                await self._execute(job)
            except Exception:  # noqa: BLE001 - a failed job must not kill the worker
                log.exception("run %s failed unexpectedly", job.run_id)
                self._set_status(job.run_id, "failed")
            finally:
                self._current = None
                self._queue.task_done()

    async def _execute(self, job: Job) -> None:
        pipeline = pipelines.get(job.pipeline_id)
        run = execute.Run(job.run_id, pipeline, execute.settings.runs_dir() / job.run_id)
        self._set_status(job.run_id, "running")

        code = await asyncio.to_thread(
            execute.execute, run, job.reads, job.reference,
            cores=job.cores, memory=job.memory, echo=False,
        )

        summary = None
        if code == 0:
            try:
                summary = await asyncio.to_thread(self._normalise, run, pipeline)
            except Exception as exc:  # noqa: BLE001 - the run itself succeeded
                log.warning("could not normalise %s: %s", job.run_id, exc)
        self._finish(job.run_id, code, summary)

    @staticmethod
    def _normalise(run: execute.Run, pipeline) -> dict:
        import yaml
        meta = yaml.safe_load((run.directory / "run.json").read_text()) or {}
        sample = (meta.get("params") or {}).get("sample_label", run.id)
        out = run.directory / "results"
        normalised = results.normalise(
            run.work, out, pipeline_id=pipeline.id,
            extra={"sample_label": sample, "run_id": run.id},
        )
        plots.render_all(normalised.per_read, out, sample=sample)
        return json.loads(normalised.summary.read_text())

    # ---------------------------------------------------------------- db helpers

    @staticmethod
    def _set_status(run_id: str, status: str, *, force: bool = False) -> None:
        """Cancelling is sticky: once a user asks for it, the container's non-zero exit must
        not be reported back as a failure. `force` is for resubmission, which is the one case
        that legitimately leaves the cancelled state."""
        with db.session() as session:
            run = session.get(db.Run, run_id)
            if run and (force or run.status != "cancelled"):
                run.status = status
                run.exit_code = None if force else run.exit_code
                session.add(run)
                session.commit()

    @staticmethod
    def _finish(run_id: str, code: int, summary: dict | None) -> None:
        with db.session() as session:
            run = session.get(db.Run, run_id)
            if not run:
                return
            was_cancelled = run.status == "cancelled"
            run.status = "cancelled" if was_cancelled else ("succeeded" if code == 0 else "failed")
            run.exit_code = code
            run.finished = datetime.now(timezone.utc)
            if summary:
                run.summary_json = json.dumps(summary)
            session.add(run)
            session.commit()


scheduler = Scheduler()
