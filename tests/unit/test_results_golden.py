"""Lock the published reference figures into CI.

These run against real outputs from the reference implementation, committed under
tests/data/golden/, so the statistics are checked on every commit in milliseconds with no
container and no 250 MB FASTQ. If a change to the normalisers moves any of these numbers, it
is a regression against Karimian et al.
"""
import csv
import json
from pathlib import Path

import pytest

from telomers import results

GOLDEN = Path(__file__).resolve().parents[1] / "data" / "golden"


@pytest.fixture
def telomere_r_work(tmp_path):
    """The reference's per-read table, laid out as pipeline A's work directory."""
    analysis = tmp_path / "r_analysis"
    analysis.mkdir()
    (analysis / "reads.table.txt").write_text((GOLDEN / "reads.table.txt").read_text())
    return tmp_path


@pytest.fixture
def telonp_work(tmp_path):
    (tmp_path / "HG002.1.NB65uq.telonp.all.csv").write_text(
        (GOLDEN / "HG002.1.NB65uq.telonp.all.csv").read_text()
    )
    return tmp_path


def test_telomere_r_reproduces_published_summary(telomere_r_work):
    rows = results.read_telomere_r(telomere_r_work)
    summary = results.summarise(rows, pipeline_id="telomere-r")
    tel = summary["telomere_bp"]
    # summary.overall_HG002.1.NB65uq.MATERNAL.csv
    assert tel["n"] == 4362
    assert tel["mean"] == 4827.2
    assert tel["median"] == 4589.5
    assert tel["sd"] == 1785.9
    assert tel["min"] == 234
    assert tel["max"] == 27258


def test_telonp_reproduces_published_summary(telonp_work):
    rows = results.read_telonp(telonp_work)
    summary = results.summarise(rows, pipeline_id="telonp")
    # HG002.1.NB65uq.telonp.summary.txt
    assert summary["reads"] == {"total": 5961, "passed": 5049, "pass_rate": 0.847}
    tel = summary["telomere_bp"]
    assert (tel["mean"], tel["median"], tel["sd"]) == (4470.6, 4378.0, 2124.6)
    assert (tel["min"], tel["max"]) == (6, 35844)
    # The sentinel breakdown, which the reference reports as err_* counts.
    assert summary["qc_breakdown"] == {"init": 674, "strand_type": 228, "fused_read": 10}


def test_per_arm_matches_the_pipelines_own_aggregation(telomere_r_work):
    """We recompute per-arm statistics from the per-read table rather than parsing the R's
    summary.by_chr, so that a filter can be applied later without re-running. That is only
    safe if the two agree."""
    rows = results.read_telomere_r(telomere_r_work)
    ours = {r["chrom_arm"]: r for r in results.per_arm(rows)}
    with (GOLDEN / "summary.by_chr.csv").open(newline="") as handle:
        theirs = {r["chrEnd"]: r for r in csv.DictReader(handle)}

    assert set(ours) == set(theirs)
    assert len(ours) == 46
    for arm, reference in theirs.items():
        for key in ("n", "median_telomere_bp", "mean_telomere_bp", "sd_telomere_bp"):
            assert abs(float(ours[arm][key]) - float(reference[key])) <= 0.05, f"{arm}.{key}"


def test_headline_biology_survives(telomere_r_work):
    """The paper's result: telomere length differs by chromosome end."""
    per_arm = {r["chrom_arm"]: r["median_telomere_bp"]
               for r in results.per_arm(results.read_telomere_r(telomere_r_work))}
    assert per_arm["chr14p"] == 2555      # shortest
    assert per_arm["chr17q"] == 2693
    assert per_arm["chr01q"] == 7124      # longest
    assert min(per_arm.values()) == per_arm["chr14p"]
    assert max(per_arm.values()) == per_arm["chr01q"]


def test_no_adapter_correction_is_applied(telomere_r_work):
    """Both pipelines measure from the read terminus with the TeloTag included, so subtracting
    it from one and not the other would introduce a bias. The value must pass through."""
    rows = results.read_telomere_r(telomere_r_work)
    with (GOLDEN / "reads.table.txt").open(newline="") as handle:
        source = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["telomere_bp"] == source[0]["Telomere.Adap"]


def test_haplotype_is_split_out_not_discarded(telomere_r_work):
    rows = results.read_telomere_r(telomere_r_work)
    assert {r["haplotype"] for r in rows} == {"MATERNAL"}
    assert all(r["chrom_arm"].startswith("chr") and not r["chrom_arm"].endswith("_M")
               for r in rows)


def test_normalise_writes_all_three_files(telonp_work, tmp_path):
    out = tmp_path / "results"
    got = results.normalise(telonp_work, out, pipeline_id="telonp")
    assert got.per_read.is_file()
    assert json.loads(got.summary.read_text())["pipeline"] == "telonp"
    # Reference-free, so there is nothing to report per arm -- the file must be absent rather
    # than empty, since the interface keys the per-arm view off its existence.
    assert got.per_arm is None


@pytest.fixture
def teloseq_work(tmp_path):
    for name in ("HG002.1.NB65uq.teloseq.perread.all.csv",
                 "HG002.1.NB65uq.teloseq.perread.csv"):
        (tmp_path / name).write_text((GOLDEN / name).read_text())
    return tmp_path


def test_teloseq_reproduces_published_summary(teloseq_work):
    rows = results.read_teloseq(teloseq_work)
    summary = results.summarise(rows, pipeline_id="teloseq")
    # HG002.1.NB65uq.teloseq.summary.txt
    assert summary["reads"]["total"] == 4234
    assert summary["reads"]["passed"] == 2521
    tel = summary["telomere_bp"]
    assert (tel["mean"], tel["median"], tel["sd"]) == (5773.7, 4665.0, 5622.5)
    assert (tel["min"], tel["max"]) == (1000, 103431)


def test_teloseq_pass_set_comes_from_the_pipeline_not_the_flags(teloseq_work):
    """Re-deriving it from is_terminal and started_correctly omits the minimum tract length,
    which reported 2,922 passing reads against the 2,521 the pipeline actually kept."""
    rows = results.read_teloseq(teloseq_work)
    flags_only = sum(1 for r in rows if r["qc_reason"] in ("", "tract_too_short"))
    assert flags_only == 2922                      # what the flags alone would have said
    assert sum(1 for r in rows if r["qc_status"] == "pass") == 2521


def test_teloseq_without_a_reference_reports_no_arms(teloseq_work):
    rows = results.read_teloseq(teloseq_work)
    assert not any(r["chrom_arm"] for r in rows)
    assert results.per_arm(rows) == []


def test_teloseq_max_tract_is_an_outlier_not_a_telomere(teloseq_work):
    """NCRF's length is end.max - start.min, so a read carrying an interstitial repeat spans
    the intervening sequence. 103 kb is not a telomere, and the metric has to be read with
    that in mind."""
    rows = results.read_teloseq(teloseq_work)
    lengths = [float(r["telomere_bp"]) for r in rows if r["qc_status"] == "pass"]
    assert max(lengths) == 103431
    assert sum(1 for v in lengths if v > 20000) / len(lengths) < 0.05
