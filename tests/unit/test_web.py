"""Web layer. Read-only, so these check routing, rendering and — mostly — that the download
routes cannot be walked out of."""
import json

import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def client(tmp_path, monkeypatch):
    from telomers import db
    from telomers.config import settings

    monkeypatch.setattr(settings, "data_root", tmp_path)
    monkeypatch.setattr(db, "_engine", None)   # engine is cached; rebind it to the temp root

    run = tmp_path / "runs" / "demo"
    (run / "results").mkdir(parents=True)
    (run / "work" / "r_analysis").mkdir(parents=True)
    (run / "run.json").write_text(json.dumps({
        "run_id": "demo", "pipeline": "telonp", "image": "telomers/pipeline-b:0.1.0",
        "params": {"sample_label": "demo"}, "inputs": {"reads": "/inputs/reads.fastq"},
    }))
    (run / "run.log").write_text("hello from the container\n")
    (run / "results" / "summary.json").write_text(json.dumps({
        "pipeline": "telonp", "sample_label": "demo",
        "reads": {"total": 10, "passed": 8, "pass_rate": 0.8},
        "qc_breakdown": {"init": 2},
        "telomere_bp": {"n": 8, "mean": 100.0, "median": 95.0, "sd": 10.0,
                        "min": 80, "max": 120, "q25": 90.0, "q75": 110.0},
    }))
    (run / "results" / "per_read.csv").write_text("read_id,telomere_bp\nr1,100\n")
    (tmp_path / "secret.txt").write_text("should never be served")

    db.reindex()
    from telomers.web.app import app
    return TestClient(app)


@pytest.mark.parametrize("path", ["/", "/runs", "/datasets", "/references", "/pipelines",
                                  "/healthz", "/runs/demo", "/runs/demo/results"])
def test_pages_render(client, path):
    assert client.get(path).status_code == 200


def test_results_page_shows_the_numbers(client):
    body = client.get("/runs/demo/results").text
    assert "95" in body and "8" in body


def test_reference_free_pipeline_explains_the_missing_per_arm_view(client):
    """An empty per-arm table would read as 'this sample has no telomeres'."""
    body = client.get("/runs/demo/results").text
    assert "never maps to" in body
    assert "<td>chr" not in body


@pytest.mark.parametrize("path", [
    "/runs/demo/work/../../secret.txt",
    "/runs/demo/results/../../secret.txt",
    "/runs/demo/work/../../../../etc/passwd",
])
def test_downloads_refuse_to_escape_the_run_directory(client, path):
    assert client.get(path).status_code in (403, 404)


def test_unknown_run_is_a_404(client):
    assert client.get("/runs/nope").status_code == 404


def test_healthz_reports_check_status(client):
    payload = client.get("/healthz").json()
    assert "ok" in payload and payload["checks"]
    assert {"name", "status", "detail"} <= set(payload["checks"][0])


# --------------------------------------------------------------------------- write path

def test_run_form_renders_fields_from_the_model(client):
    body = client.get("/runs/new?pipeline=telomere-r").text
    assert 'name="barcode_name"' in body and 'name="sample_label"' in body
    # chr_arm_ln comes from the reference, so it must not be an editable input.
    assert 'name="chr_arm_ln"' not in body
    assert "from the selected reference" in body


def test_reference_free_pipeline_form_has_no_reference_selector(client):
    assert 'name="reference_id"' not in client.get("/runs/new?pipeline=telonp").text


def test_runs_new_is_not_mistaken_for_a_run_id(client):
    """/runs/{run_id} is also registered; order matters."""
    assert client.get("/runs/new").status_code == 200


def test_path_registration_is_disabled_until_an_allowlist_is_set(client):
    """Registering by path reads an arbitrary local file, so it defaults to off rather than
    to permissive."""
    response = client.post("/datasets", data={"path": "/etc/passwd"}, follow_redirects=False)
    assert response.status_code == 303
    assert "disabled" in response.headers["location"]


def test_registering_a_path_outside_the_allowlist_is_refused(client, tmp_path, monkeypatch):
    from telomers.config import settings
    monkeypatch.setattr(settings, "import_dirs", [tmp_path / "allowed"])
    response = client.post("/datasets", data={"path": "/etc/passwd"}, follow_redirects=False)
    assert "not+under+any+allowed" in response.headers["location"]


def test_submitting_without_a_dataset_is_refused(client):
    response = client.post("/runs", data={"pipeline": "telonp"}, follow_redirects=False)
    assert "choose+a+dataset" in response.headers["location"]


def test_mapping_pipeline_refuses_to_start_without_a_reference(client, tmp_path):
    from telomers import db
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    with db.session() as session:
        session.add(db.Dataset(id="d1", name="d", path=str(reads), sha256="x", size_bytes=1))
        session.commit()
    response = client.post("/runs", data={"pipeline": "telomere-r", "dataset_id": "d1"},
                           follow_redirects=False)
    assert "needs+a+reference" in response.headers["location"]


def test_status_fragment_polls_only_while_live(client):
    from telomers import db
    # The demo run is finished, so the fragment must not carry hx-trigger -- that is how
    # polling stops without any client-side logic.
    assert "hx-trigger" not in client.get("/runs/demo/status").text

    with db.session() as session:
        run = session.get(db.Run, "demo")
        run.status = "running"
        session.add(run)
        session.commit()
    assert "hx-trigger" in client.get("/runs/demo/status").text
