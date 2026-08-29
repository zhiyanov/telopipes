"""Environment checks.

The point is to fail with something actionable before a run starts, rather than deep inside a
container. The x86-emulation check earns its place: this project originally planned to run
amd64 images under FEX, and that turned out to be impossible on a 16 KB-page kernel -- see
docs/ENVIRONMENT.md.
"""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
from dataclasses import dataclass

from .config import settings
from .pipelines import registry

OK, WARN, FAIL = "ok", "warn", "fail"


@dataclass
class Check:
    name: str
    status: str
    detail: str


def _run(*cmd: str, timeout: int = 30) -> tuple[int, str]:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return proc.returncode, (proc.stdout + proc.stderr).strip()
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 1, str(exc)


def run_checks() -> list[Check]:
    checks: list[Check] = []

    engine = settings.container_engine
    path = shutil.which(engine)
    if path is None:
        checks.append(Check("container engine", FAIL, f"{engine} not found on PATH"))
        return checks
    code, out = _run(engine, "--version")
    checks.append(Check("container engine", OK if code == 0 else FAIL, out.splitlines()[0] if out else path))

    checks.append(Check("host architecture", OK, f"{platform.machine()}, page size {os.sysconf('SC_PAGE_SIZE')}"))

    # A 16 KB page size makes x86_64 emulation impossible, so a cross-platform setting here is
    # a configuration error worth catching early rather than a slow path.
    if settings.platform:
        page = os.sysconf("SC_PAGE_SIZE")
        want_x86 = "amd64" in settings.platform or "x86" in settings.platform
        if want_x86 and platform.machine() not in ("x86_64", "amd64") and page != 4096:
            checks.append(Check(
                "platform override", FAIL,
                f"TELOMERS_PLATFORM={settings.platform} requests x86 emulation, but this host "
                f"has a {page}-byte page size. x86_64 requires 4 KB pages; FEX and qemu-user "
                "both fail. See docs/ENVIRONMENT.md.",
            ))
        else:
            checks.append(Check("platform override", OK, settings.platform))

    code, out = _run("getenforce")
    if code == 0:
        selinux = out.strip()
        checks.append(Check(
            "SELinux", OK,
            f"{selinux} (mounts are labelled ,z)" if selinux == "Enforcing" else selinux,
        ))

    try:
        reg = registry()
    except Exception as exc:  # noqa: BLE001 - surfaced to the user verbatim
        checks.append(Check("pipelines", FAIL, str(exc)))
        return checks
    checks.append(Check("pipelines", OK, ", ".join(sorted(reg))))

    code, out = _run(engine, "images", "--format", "{{.Repository}}:{{.Tag}}")
    available = set(out.split()) if code == 0 else set()
    for pid, pipeline in sorted(reg.items()):
        present = any(image.endswith(pipeline.image) for image in available)
        checks.append(Check(
            f"image {pid}", OK if present else WARN,
            pipeline.image if present else f"{pipeline.image} not built -- run: telomers images build",
        ))

    root = settings.data_root
    try:
        root.mkdir(parents=True, exist_ok=True)
        usage = shutil.disk_usage(root)
        free_gb = usage.free / 1024**3
        checks.append(Check(
            "data root", OK if free_gb > 10 else WARN,
            f"{root} ({free_gb:.0f} GB free)",
        ))
    except OSError as exc:
        checks.append(Check("data root", FAIL, f"{root}: {exc}"))

    return checks
