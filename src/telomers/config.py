"""Application settings.

Everything is overridable through TELOMERS_* environment variables or a .env file. The whole
configuration surface is deliberately small: this runs on one machine for a handful of people.
"""
from __future__ import annotations

from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

_REPO_ROOT = Path(__file__).resolve().parents[2]


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="TELOMERS_", env_file=".env", extra="ignore"
    )

    #: Where datasets, references and runs live.
    data_root: Path = _REPO_ROOT / "data"
    #: Directory holding one subdirectory per pipeline (Snakefile, Dockerfile, meta.toml).
    pipelines_dir: Path = _REPO_ROOT / "pipelines"

    container_engine: str = "podman"
    #: Empty means "whatever the host is". Set explicitly only to cross-build; note that x86
    #: emulation does not work on a 16 KB-page kernel -- see docs/ENVIRONMENT.md.
    platform: str = ""

    #: Defaults for a run; per-run values override these.
    threads: int = 4
    memory: str = "10g"

    #: Directories a dataset may be registered from in place. Empty disables path import,
    #: which is otherwise an arbitrary local-file-read primitive.
    import_dirs: list[Path] = []

    def runs_dir(self) -> Path:
        return self.data_root / "runs"


settings = Settings()
