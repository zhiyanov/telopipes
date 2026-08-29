"""On-disk layout.

The only module that knows where things live. Every object is a directory holding its payload
plus a meta.json, which makes the SQLite index rebuildable by walking the data root -- see
db.reindex(). That property is why there are no migrations.

    $TELOMERS_DATA_ROOT/
      datasets/<id>/    reads.fastq[.gz]   meta.json
      references/<id>/  reference.fasta    meta.json
      runs/<id>/        run.json  run.log  work/  results/
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from .config import settings

_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def slugify(name: str, maxlen: int = 48) -> str:
    return _SAFE.sub("-", name).strip("-")[:maxlen] or "unnamed"


def unique_id(prefix: str, parent: Path) -> str:
    """A readable, collision-free directory name."""
    stamp = f"{datetime.now(timezone.utc):%Y%m%d-%H%M%S}"
    candidate = f"{stamp}-{prefix}"
    n = 1
    while (parent / candidate).exists():
        n += 1
        candidate = f"{stamp}-{prefix}-{n}"
    return candidate


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def datasets_dir() -> Path:
    return settings.data_root / "datasets"


def references_dir() -> Path:
    return settings.data_root / "references"


def runs_dir() -> Path:
    return settings.data_root / "runs"


def write_meta(directory: Path, meta: dict) -> None:
    (directory / "meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")


def read_meta(directory: Path) -> dict | None:
    path = directory / "meta.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return None


@dataclass(frozen=True)
class ImportResult:
    path: Path
    sha256: str
    size: int
    linked: bool


def ingest(source: Path, destination: Path, *, allow_link: bool = True) -> ImportResult:
    """Bring a file into the data root, hardlinking when possible.

    A hardlink means registering the 3 GB of FASTQs already on this machine costs nothing.
    It also means the file is shared, not copied -- acceptable because everything under the
    data root is treated as read-only once written.
    """
    source = source.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    linked = False
    if allow_link:
        try:
            os.link(source, destination)
            linked = True
        except OSError:
            pass  # different filesystem, or the link limit -- fall back to copying
    if not linked:
        shutil.copy2(source, destination)
    return ImportResult(destination, sha256_of(destination), destination.stat().st_size, linked)


def check_import_allowed(source: Path) -> None:
    """Registering a path is an arbitrary local-file-read primitive, so it is gated.

    Empty import_dirs disables it entirely rather than defaulting to permissive.
    """
    source = source.resolve()
    allowed = [Path(d).resolve() for d in settings.import_dirs]
    if not allowed:
        raise PermissionError(
            "registering files by path is disabled. Set TELOMERS_IMPORT_DIRS to the "
            "directories that may be read, e.g. TELOMERS_IMPORT_DIRS='[\"/data/nanopore\"]'"
        )
    if not any(source == d or d in source.parents for d in allowed):
        raise PermissionError(
            f"{source} is not under any allowed import directory: "
            + ", ".join(str(d) for d in allowed)
        )


def register_dataset(source: Path, *, name: str | None = None, link: bool = True,
                     count_reads: bool = True) -> dict:
    """Bring a FASTQ into the data root and describe it."""
    from . import validate

    source = Path(source)
    if not source.is_file():
        raise FileNotFoundError(source)

    info = validate.sniff_fastq(source, limit=None if count_reads else 1)
    if not info.ok:
        raise ValueError(f"{source.name}: " + "; ".join(info.errors))

    label = name or source.name
    directory = datasets_dir() / unique_id(slugify(Path(label).stem), datasets_dir())
    suffix = ".fastq.gz" if source.name.endswith(".gz") else ".fastq"
    result = ingest(source, directory / f"reads{suffix}", allow_link=link)

    meta = {
        "id": directory.name,
        "name": label,
        "path": str(result.path),
        "sha256": result.sha256,
        "size_bytes": result.size,
        "n_reads": info.n_reads if count_reads else None,
        "total_bases": info.total_bases if count_reads else None,
        "guessed_barcode": validate.guess_barcode(source.name),
        "source": str(source.resolve()),
        "linked": result.linked,
        "created": datetime.now(timezone.utc).isoformat(),
    }
    write_meta(directory, meta)
    return meta


def register_reference(source: Path, *, name: str | None = None, link: bool = True) -> dict:
    """Validate a cut reference and bring it into the data root.

    Rejected outright if it does not validate: an unusable reference caught here is a clear
    error message, whereas the same reference caught nowhere produces plausible-looking but
    wrong q-arm telomere lengths. See validate.validate_reference.
    """
    from . import validate

    source = Path(source)
    if not source.is_file():
        raise FileNotFoundError(source)

    info = validate.validate_reference(source)
    if not info.ok:
        raise ValueError(f"{source.name} is not a usable cut reference:\n  - "
                         + "\n  - ".join(info.errors))

    label = name or source.name
    directory = references_dir() / unique_id(slugify(Path(label).stem), references_dir())
    suffix = ".fasta.gz" if source.name.endswith(".gz") else ".fasta"
    result = ingest(source, directory / f"reference{suffix}", allow_link=link)

    meta = {
        "id": directory.name,
        "name": label,
        "path": str(result.path),
        "sha256": result.sha256,
        "size_bytes": result.size,
        "n_records": info.n_records,
        "arm_length_bp": info.arm_length_bp,
        "haplotypes": ",".join(info.haplotypes),
        "warnings": info.warnings,
        "source": str(source.resolve()),
        "linked": result.linked,
        "created": datetime.now(timezone.utc).isoformat(),
    }
    write_meta(directory, meta)
    return meta
