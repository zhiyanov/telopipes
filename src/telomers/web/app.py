"""FastAPI application.

Read-only for now: it browses whatever the CLI has produced. Submitting runs comes next.

Server-rendered Jinja with no client-side framework and no build step. For a lab tool whose
pages are tables, figures and download links, a SPA would be all cost and no benefit -- and it
would put a node toolchain between the biologists and their results.
"""
from __future__ import annotations

import json
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, HTMLResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from sqlmodel import select

from .. import db, doctor, pipelines, storage
from ..config import settings

HERE = Path(__file__).resolve().parent
templates = Jinja2Templates(directory=str(HERE / "templates"))

app = FastAPI(title="telomers", docs_url=None, redoc_url=None)
app.mount("/static", StaticFiles(directory=str(HERE / "static")), name="static")


# --------------------------------------------------------------------------- helpers

def _fmt_bytes(n: int | None) -> str:
    if not n:
        return "—"
    value = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:,.0f} {unit}" if unit == "B" else f"{value:,.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TB"


def _fmt_int(n) -> str:
    return f"{int(n):,}" if n not in (None, "") else "—"


templates.env.filters["bytes"] = _fmt_bytes
templates.env.filters["count"] = _fmt_int


def _run_dir(run_id: str) -> Path:
    """Resolve a run directory, refusing anything that escapes the runs root."""
    root = storage.runs_dir().resolve()
    directory = (root / run_id).resolve()
    if not (directory == root or root in directory.parents) or not directory.is_dir():
        raise HTTPException(404, f"no such run: {run_id}")
    return directory


def _safe_child(base: Path, relative: str) -> Path:
    """Resolve a path inside `base`, rejecting traversal.

    Every download route funnels through here; the check is on the *resolved* path, so
    symlinks cannot be used to step outside either.
    """
    base = base.resolve()
    target = (base / relative).resolve()
    if base not in target.parents:
        raise HTTPException(403, "path outside the run directory")
    if not target.is_file():
        raise HTTPException(404, relative)
    return target


def _load_summary(directory: Path) -> dict | None:
    path = directory / "results" / "summary.json"
    return json.loads(path.read_text()) if path.is_file() else None


def _read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    import csv
    if not path.is_file():
        return [], []
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    return (rows[0], rows[1:]) if rows else ([], [])


# --------------------------------------------------------------------------- pages

@app.get("/", response_class=HTMLResponse)
def dashboard(request: Request):
    with db.session() as session:
        runs = session.exec(select(db.Run).order_by(db.Run.created.desc())).all()
        n_datasets = len(session.exec(select(db.Dataset)).all())
        n_references = len(session.exec(select(db.Reference)).all())
    import shutil
    settings.data_root.mkdir(parents=True, exist_ok=True)
    usage = shutil.disk_usage(settings.data_root)
    return templates.TemplateResponse(request, "index.html", {
        "runs": runs[:8],
        "n_runs": len(runs),
        "n_datasets": n_datasets,
        "n_references": n_references,
        "free_bytes": usage.free,
        "pipelines": sorted(pipelines.registry().values(), key=lambda p: p.id),
    })


@app.get("/runs", response_class=HTMLResponse)
def runs_list(request: Request):
    with db.session() as session:
        runs = session.exec(select(db.Run).order_by(db.Run.created.desc())).all()
    return templates.TemplateResponse(request, "runs.html", {"runs": runs})


@app.get("/runs/{run_id}", response_class=HTMLResponse)
def run_detail(request: Request, run_id: str):
    directory = _run_dir(run_id)
    with db.session() as session:
        run = session.get(db.Run, run_id)

    meta = {}
    run_json = directory / "run.json"
    if run_json.is_file():
        import yaml
        meta = yaml.safe_load(run_json.read_text()) or {}

    log = ""
    if (directory / "run.log").is_file():
        log = (directory / "run.log").read_text()[-8000:]

    work = directory / "work"
    pipeline_figures = sorted(p.relative_to(work).as_posix()
                              for p in work.rglob("*.png")) if work.is_dir() else []
    return templates.TemplateResponse(request, "run_detail.html", {
        "run_id": run_id,
        "run": run,
        "meta": meta,
        "summary": _load_summary(directory),
        "log": log,
        "n_figures": len(pipeline_figures),
    })


@app.get("/runs/{run_id}/results", response_class=HTMLResponse)
def run_results(request: Request, run_id: str):
    directory = _run_dir(run_id)
    summary = _load_summary(directory)
    if summary is None:
        raise HTTPException(404, "this run has no canonical results yet")

    results = directory / "results"
    # Deliberate order: the bulk distribution is the headline, the per-arm breakdown is the
    # reason to use pipeline A at all, and the funnel is diagnostic. Alphabetical would put
    # the diagnostic in the middle.
    order = ["distribution.png", "by_arm.png", "read_funnel.png"]
    present = {p.name for p in results.glob("*.png")}
    figures = [n for n in order if n in present] + sorted(present - set(order))
    arm_header, arm_rows = _read_csv(results / "per_arm.csv")

    work = directory / "work"
    pipeline_figures = sorted(p.relative_to(work).as_posix()
                              for p in work.rglob("*.png")) if work.is_dir() else []

    downloads = [(p.name, p.stat().st_size)
                 for p in sorted(results.iterdir()) if p.suffix in (".csv", ".json")]
    return templates.TemplateResponse(request, "results.html", {
        "run_id": run_id,
        "summary": summary,
        "figures": figures,
        "arm_header": arm_header,
        "arm_rows": arm_rows,
        "pipeline_figures": pipeline_figures,
        "downloads": downloads,
    })


@app.get("/datasets", response_class=HTMLResponse)
def datasets_list(request: Request):
    with db.session() as session:
        rows = session.exec(select(db.Dataset).order_by(db.Dataset.created.desc())).all()
    return templates.TemplateResponse(request, "datasets.html", {"datasets": rows})


@app.get("/references", response_class=HTMLResponse)
def references_list(request: Request):
    with db.session() as session:
        rows = session.exec(select(db.Reference).order_by(db.Reference.created.desc())).all()
    return templates.TemplateResponse(request, "references.html", {"references": rows})


@app.get("/pipelines", response_class=HTMLResponse)
def pipelines_list(request: Request):
    return templates.TemplateResponse(request, "pipelines.html", {
        "pipelines": sorted(pipelines.registry().values(), key=lambda p: p.id),
    })


# --------------------------------------------------------------------------- files

@app.get("/runs/{run_id}/results/{name}")
def results_file(run_id: str, name: str):
    target = _safe_child(_run_dir(run_id) / "results", name)
    inline = target.suffix in (".png", ".json")
    return FileResponse(target, filename=None if inline else target.name)


@app.get("/runs/{run_id}/work/{path:path}")
def work_file(run_id: str, path: str):
    target = _safe_child(_run_dir(run_id) / "work", path)
    return FileResponse(target, filename=None if target.suffix == ".png" else target.name)


@app.get("/runs/{run_id}/log", response_class=PlainTextResponse)
def run_log(run_id: str):
    path = _run_dir(run_id) / "run.log"
    if not path.is_file():
        raise HTTPException(404, "no log")
    return path.read_text()


@app.get("/healthz")
def healthz():
    checks = doctor.run_checks()
    return {
        "ok": all(c.status != doctor.FAIL for c in checks),
        "checks": [{"name": c.name, "status": c.status, "detail": c.detail} for c in checks],
    }
