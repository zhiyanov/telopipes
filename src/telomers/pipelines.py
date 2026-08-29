"""Pipeline registry.

A pipeline is a directory under `pipelines/` holding a Snakefile, a Dockerfile and a meta.toml.
meta.toml carries display metadata only -- parameters live in params.py -- so adding a pipeline
is a directory plus one entry in PARAMS_MODELS.
"""
from __future__ import annotations

import tomllib
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from pydantic import BaseModel

from .config import settings
from .params import PARAMS_MODELS


@dataclass(frozen=True)
class Pipeline:
    id: str
    version: str
    title: str
    summary: str
    citation: str
    upstream: str
    image: str
    requires_reference: bool
    provides: tuple[str, ...]
    directory: Path

    @property
    def params_model(self) -> type[BaseModel]:
        return PARAMS_MODELS[self.id]

    @property
    def snakefile(self) -> Path:
        return self.directory / "Snakefile"

    def has(self, capability: str) -> bool:
        """Whether this pipeline produces a given kind of result.

        Drives the UI: pipeline B has no per_arm output, and that tab should be absent rather
        than empty -- an empty per-arm chart reads as "this sample has no telomeres".
        """
        return capability in self.provides


def _load(directory: Path) -> Pipeline:
    meta = tomllib.loads((directory / "meta.toml").read_text())["pipeline"]
    missing = {"id", "title", "image", "requires_reference", "provides"} - set(meta)
    if missing:
        raise ValueError(f"{directory/'meta.toml'} is missing: {', '.join(sorted(missing))}")
    if meta["id"] not in PARAMS_MODELS:
        raise ValueError(
            f"pipeline {meta['id']!r} has no parameter model; add one to params.PARAMS_MODELS"
        )
    if not (directory / "Snakefile").is_file():
        raise ValueError(f"{directory} has no Snakefile")
    return Pipeline(
        id=meta["id"],
        version=meta.get("version", "0"),
        title=meta["title"],
        summary=meta.get("summary", "").strip(),
        citation=meta.get("citation", ""),
        upstream=meta.get("upstream", ""),
        image=meta["image"],
        requires_reference=bool(meta["requires_reference"]),
        provides=tuple(meta["provides"]),
        directory=directory,
    )


@cache
def registry() -> dict[str, Pipeline]:
    """All discoverable pipelines, keyed by id. Validated eagerly so a malformed one fails
    at startup rather than when someone tries to run it."""
    root = settings.pipelines_dir
    if not root.is_dir():
        raise FileNotFoundError(f"pipelines directory not found: {root}")
    found = {}
    for directory in sorted(p for p in root.iterdir() if (p / "meta.toml").is_file()):
        pipeline = _load(directory)
        found[pipeline.id] = pipeline
    if not found:
        raise FileNotFoundError(f"no pipelines found under {root}")
    return found


def get(pipeline_id: str) -> Pipeline:
    reg = registry()
    if pipeline_id not in reg:
        raise KeyError(f"unknown pipeline {pipeline_id!r}; available: {', '.join(sorted(reg))}")
    return reg[pipeline_id]
