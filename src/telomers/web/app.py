"""FastAPI application.

Read-only for now: it browses whatever the CLI has produced. Submitting runs comes next.

Server-rendered Jinja with no client-side framework and no build step. For a lab tool whose
pages are tables, figures and download links, a SPA would be all cost and no benefit -- and it
would put a node toolchain between the biologists and their results.
"""
from __future__ import annotations

import json
from pathlib import Path

from contextlib import asynccontextmanager

from fastapi import FastAPI, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, PlainTextResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import ValidationError
from sqlmodel import select

from .. import db, doctor, execute, forms, jobs, pipelines, storage, validate
from ..config import settings

HERE = Path(__file__).resolve().parent
templates = Jinja2Templates(directory=str(HERE / "templates"))

@asynccontextmanager
async def lifespan(_: FastAPI):
    # Any run still marked running when the process died cannot be adopted, so it is marked
    # interrupted and offered a Resume. Snakemake state lives in the run directory, so
    # resuming re-executes only what did not finish.
    await jobs.scheduler.start()
    yield
    await jobs.scheduler.stop()


app = FastAPI(title="telomers", docs_url=None, redoc_url=None, lifespan=lifespan)
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


@app.get("/runs/new", response_class=HTMLResponse)
def run_new(request: Request, pipeline: str = ""):
    registry = pipelines.registry()
    chosen = registry.get(pipeline)
    with db.session() as session:
        datasets = session.exec(select(db.Dataset).order_by(db.Dataset.created.desc())).all()
        references = session.exec(select(db.Reference).order_by(db.Reference.created.desc())).all()
    return templates.TemplateResponse(request, "run_new.html", {
        "pipelines": sorted(registry.values(), key=lambda p: p.id),
        "chosen": chosen,
        "fields": forms.describe(chosen.params_model) if chosen else [],
        "datasets": datasets,
        "references": references,
    })


# Registered before /runs/{run_id}: FastAPI matches in definition order, so putting
# this later would make "new" look like a run id.
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


# --------------------------------------------------------------------------- write path

def _flash(request: Request, message: str, kind: str = "ok") -> RedirectResponse:
    """Carry a one-line result across a redirect via the query string.

    Simpler than session middleware for a single-user lab tool, and it survives a refresh.
    """
    from urllib.parse import urlencode
    target = request.headers.get("referer", "/")
    joiner = "&" if "?" in target else "?"
    return RedirectResponse(f"{target}{joiner}{urlencode({'msg': message, 'kind': kind})}",
                            status_code=303)


@app.post("/datasets")
async def datasets_create(request: Request,
                          path: str = Form(""),
                          name: str = Form(""),
                          upload: UploadFile | None = None):
    try:
        if path:
            storage.check_import_allowed(Path(path))
            meta = storage.register_dataset(Path(path), name=name or None)
        elif upload is not None and upload.filename:
            staged = _stage_upload(upload)
            try:
                meta = storage.register_dataset(staged, name=name or upload.filename, link=False)
            finally:
                staged.unlink(missing_ok=True)
        else:
            return _flash(request, "give a server path or choose a file", "bad")
    except (PermissionError, ValueError, FileNotFoundError) as exc:
        return _flash(request, str(exc), "bad")

    with db.session() as session:
        session.add(db.from_meta(db.Dataset, meta))
        session.commit()
    return _flash(request, f"registered {meta['name']} ({meta['n_reads']:,} reads)")


@app.post("/references")
async def references_create(request: Request,
                            path: str = Form(""),
                            name: str = Form(""),
                            cut: str = Form(""),
                            haplotype: str = Form("all"),
                            upload: UploadFile | None = None):
    opts = {"cut": bool(cut), "haplotype": haplotype if haplotype in
            ("all", "MATERNAL", "PATERNAL") else "all"}
    try:
        if path:
            storage.check_import_allowed(Path(path))
            meta = storage.register_reference(Path(path), name=name or None, **opts)
        elif upload is not None and upload.filename:
            staged = _stage_upload(upload)
            try:
                meta = storage.register_reference(staged, name=name or upload.filename,
                                                  link=False, **opts)
            finally:
                staged.unlink(missing_ok=True)
        else:
            return _flash(request, "give a server path or choose a file", "bad")
    except (PermissionError, ValueError, FileNotFoundError, RuntimeError) as exc:
        return _flash(request, str(exc), "bad")

    with db.session() as session:
        session.add(db.from_meta(db.Reference, meta))
        session.commit()
    return _flash(request, f"registered {meta['name']} "
                           f"({meta['n_records']} records, arm {meta['arm_length_bp']:,} bp)")


def _stage_upload(upload: UploadFile) -> Path:
    """Write an uploaded file to a staging path.

    Starlette spools to a temporary file on disk past ~1 MB, so this does not hold a 250 MB
    FASTQ in memory -- but it does mean the bytes are written twice. Registering by path
    avoids that entirely, and is the route that will actually be used for data already on
    this machine.
    """
    import shutil as _shutil
    staging = settings.data_root / "staging"
    staging.mkdir(parents=True, exist_ok=True)
    target = staging / Path(upload.filename).name
    with target.open("wb") as handle:
        _shutil.copyfileobj(upload.file, handle, length=1 << 20)
    return target


@app.post("/runs")
async def run_create(request: Request):
    form = dict(await request.form())
    pipeline_id = str(form.get("pipeline", ""))
    try:
        pipeline = pipelines.get(pipeline_id)
    except KeyError as exc:
        raise HTTPException(400, str(exc)) from exc

    with db.session() as session:
        dataset = session.get(db.Dataset, str(form.get("dataset_id", "")))
        reference = session.get(db.Reference, str(form.get("reference_id", "")))
    if dataset is None:
        return _flash(request, "choose a dataset", "bad")
    if pipeline.requires_reference and reference is None:
        return _flash(request, f"{pipeline.title} needs a reference genome", "bad")

    values = forms.coerce(pipeline.params_model, form)
    values.setdefault("sample_label", Path(dataset.name).stem.replace(" ", "_"))
    # The reference determines this; a value that disagrees with the file would silently
    # corrupt q-arm lengths, so the form never lets it be typed and we set it here.
    if reference is not None and "chr_arm_ln" in pipeline.params_model.model_fields:
        values["chr_arm_ln"] = reference.arm_length_bp
    try:
        params = pipeline.params_model(**values)
    except ValidationError as exc:
        first = exc.errors()[0]
        where = ".".join(str(p) for p in first["loc"]) or "parameters"
        return _flash(request, f"{where}: {first['msg']}", "bad")

    reads = Path(dataset.path)
    reference_path = Path(reference.path) if reference else None
    run = execute.prepare(pipeline, params, reads, reference_path)

    with db.session() as session:
        session.add(db.Run(
            id=run.id, pipeline=pipeline.id, pipeline_version=pipeline.version,
            status="created", dataset_id=dataset.id,
            reference_id=reference.id if reference else None,
            params_json=json.dumps(params.model_dump()),
        ))
        session.commit()

    jobs.scheduler.submit(jobs.Job(run.id, pipeline.id, reads, reference_path))
    return RedirectResponse(f"/runs/{run.id}", status_code=303)


@app.post("/runs/{run_id}/cancel")
async def run_cancel(run_id: str):
    _run_dir(run_id)
    await jobs.scheduler.cancel(run_id)
    return RedirectResponse(f"/runs/{run_id}", status_code=303)


@app.post("/runs/{run_id}/resume")
async def run_resume(run_id: str):
    directory = _run_dir(run_id)
    with db.session() as session:
        run = session.get(db.Run, run_id)
        if run is None:
            raise HTTPException(404, run_id)
        dataset = session.get(db.Dataset, run.dataset_id) if run.dataset_id else None
        reference = session.get(db.Reference, run.reference_id) if run.reference_id else None
        pipeline_id = run.pipeline

    import yaml
    meta = yaml.safe_load((directory / "run.json").read_text()) or {}
    reads = Path(dataset.path) if dataset else Path(meta["inputs"]["reads"])
    ref = Path(reference.path) if reference else (
        Path(meta["inputs"]["reference"]) if meta.get("inputs", {}).get("reference") else None
    )
    jobs.scheduler.submit(jobs.Job(run_id, pipeline_id, reads, ref))
    return RedirectResponse(f"/runs/{run_id}", status_code=303)


@app.post("/runs/{run_id}/delete")
async def run_delete(run_id: str):
    import shutil as _shutil
    directory = _run_dir(run_id)
    await jobs.scheduler.cancel(run_id)
    _shutil.rmtree(directory, ignore_errors=True)
    with db.session() as session:
        run = session.get(db.Run, run_id)
        if run:
            session.delete(run)
            session.commit()
    return RedirectResponse("/runs", status_code=303)


@app.get("/runs/{run_id}/status", response_class=HTMLResponse)
def run_status(request: Request, run_id: str):
    """HTMX fragment, polled while a run is live.

    The server stops emitting the polling attribute once the run reaches a terminal state, so
    polling stops on its own with no client-side logic.
    """
    directory = _run_dir(run_id)
    with db.session() as session:
        run = session.get(db.Run, run_id)
    log_path = directory / "run.log"
    return templates.TemplateResponse(request, "_status.html", {
        "run_id": run_id,
        "run": run,
        "summary": _load_summary(directory),
        "log": log_path.read_text()[-8000:] if log_path.is_file() else "",
        "live": run is not None and run.status in ("queued", "running"),
        "is_current": jobs.scheduler.current == run_id,
    })
