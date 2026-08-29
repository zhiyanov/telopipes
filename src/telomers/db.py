"""SQLite index over the data root.

The database is a *cache*, not the record. Every fact also lives on disk in the meta.json or
run.json next to the object it describes, so `reindex()` can rebuild it by walking the data
root. That is what lets this project skip migrations entirely without risking run history: if
the schema changes, delete the file and reindex.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import Engine
from sqlmodel import Field, Session, SQLModel, create_engine, select

__all__ = ["Dataset", "Reference", "Run", "engine", "session", "reindex", "select",
           "from_meta"]

from . import storage
from .config import settings


def _from_meta(model: type[SQLModel], meta: dict) -> SQLModel:
    """Build a row from an on-disk meta.json.

    meta.json stores timestamps as ISO strings, which is what makes it readable and portable;
    the models want datetimes. Coercing in one place keeps every call site -- registration and
    reindex alike -- from having to remember.
    """
    fields = model.model_fields
    values = {}
    for key, value in meta.items():
        if key not in fields:
            continue
        annotation = fields[key].annotation
        if isinstance(value, str) and annotation in (datetime, datetime | None):
            try:
                value = datetime.fromisoformat(value)
            except ValueError:
                continue
        values[key] = value
    return model(**values)


class Dataset(SQLModel, table=True):
    id: str = Field(primary_key=True)
    name: str
    path: str
    sha256: str
    size_bytes: int
    n_reads: int | None = None
    total_bases: int | None = None
    guessed_barcode: str | None = None
    created: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))


class Reference(SQLModel, table=True):
    id: str = Field(primary_key=True)
    name: str
    path: str
    sha256: str
    size_bytes: int
    n_records: int
    #: Derived from the file, never supplied by a user -- see validate.validate_reference.
    arm_length_bp: int
    haplotypes: str = ""
    created: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))


class Run(SQLModel, table=True):
    id: str = Field(primary_key=True)
    pipeline: str
    pipeline_version: str = ""
    status: str = "created"          # created | running | succeeded | failed | interrupted
    dataset_id: str | None = None
    reference_id: str | None = None
    params_json: str = "{}"
    summary_json: str | None = None
    exit_code: int | None = None
    created: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    finished: datetime | None = None

    @property
    def params(self) -> dict:
        return json.loads(self.params_json)

    @property
    def summary(self) -> dict | None:
        return json.loads(self.summary_json) if self.summary_json else None


_engine: Engine | None = None


def engine() -> Engine:
    global _engine
    if _engine is None:
        settings.data_root.mkdir(parents=True, exist_ok=True)
        path = settings.data_root / "telomers.sqlite3"
        _engine = create_engine(f"sqlite:///{path}", echo=False)
        SQLModel.metadata.create_all(_engine)
    return _engine


def session() -> Session:
    return Session(engine())


def reindex() -> dict[str, int]:
    """Rebuild the index from what is on disk. Safe to run at any time."""
    counts = {"datasets": 0, "references": 0, "runs": 0}
    with session() as db:
        for row in db.exec(select(Dataset)).all():
            db.delete(row)
        for row in db.exec(select(Reference)).all():
            db.delete(row)
        for row in db.exec(select(Run)).all():
            db.delete(row)
        db.commit()

        for directory in sorted(_subdirs(storage.datasets_dir())):
            meta = storage.read_meta(directory)
            if meta:
                db.add(_from_meta(Dataset, meta))
                counts["datasets"] += 1

        for directory in sorted(_subdirs(storage.references_dir())):
            meta = storage.read_meta(directory)
            if meta:
                db.add(_from_meta(Reference, meta))
                counts["references"] += 1

        for directory in sorted(_subdirs(storage.runs_dir())):
            run = _run_from_disk(directory)
            if run:
                db.add(run)
                counts["runs"] += 1
        db.commit()
    return counts


def _subdirs(parent: Path) -> list[Path]:
    return [p for p in parent.iterdir() if p.is_dir()] if parent.is_dir() else []


def _run_from_disk(directory: Path) -> Run | None:
    path = directory / "run.json"
    if not path.is_file():
        return None
    try:
        import yaml
        meta = yaml.safe_load(path.read_text())
    except Exception:  # noqa: BLE001 - a corrupt run.json should not abort a reindex
        return None
    if not isinstance(meta, dict) or "run_id" not in meta:
        return None

    summary_path = directory / "results" / "summary.json"
    # A run whose outputs exist has clearly finished; without a status file that is the best
    # signal available, and it is why reindex can recover from a lost database.
    status = "succeeded" if summary_path.is_file() else "interrupted"
    return Run(
        id=meta["run_id"],
        pipeline=meta.get("pipeline", ""),
        pipeline_version=meta.get("pipeline_version", ""),
        status=status,
        params_json=json.dumps(meta.get("params", {})),
        summary_json=summary_path.read_text() if summary_path.is_file() else None,
    )


#: Public alias -- registration in cli.py needs the same coercion as reindex.
from_meta = _from_meta
