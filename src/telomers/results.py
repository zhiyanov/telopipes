"""Normalise each pipeline's output into one canonical form.

Both pipelines are reduced to a per-read table with the same columns, so the interface can be
written once. Only stdlib is used -- a few thousand rows does not justify pandas, and keeping
it out means the host dependency list stays short.

**No adapter correction is applied, in either pipeline.** It is tempting to assume pipeline A's
`Telomere+Adap` needs the ligated TeloTag subtracted before it can be compared with pipeline
B's length, but both already measure the same quantity: bases from the read terminus to the
telomere boundary, tag included. Pipeline A's is `softClip - mapStart + 1`, anchored at the
terminus. Pipeline B's `getTeloNPBoundary` returns `boundaryPoint = x * windowStep` from a
window loop that starts at index 0 (C strand) or at the last base (G strand), so the tag falls
inside the returned index too -- it is simply not telomere-like, and sits outside the window
that decides where the boundary is. Subtracting 24 bp from one and not the other would
*introduce* a bias. See docs/PLAN.md section 5.
"""
from __future__ import annotations

import csv
import json
import re
import statistics
from dataclasses import dataclass
from pathlib import Path

#: The canonical per-read schema. Pipeline B leaves the mapping columns empty.
PER_READ_COLUMNS = [
    "read_id",
    "telomere_bp",
    "read_length_bp",
    "chrom_arm",
    "haplotype",
    "strand",
    "mapq",
    "qc_status",
    "qc_reason",
]

#: TeloBP's failure sentinels, returned in place of a length.
TELONP_SENTINELS = {
    -1: "init",
    -10: "fused_read",
    -20: "strand_type",
    -1000: "seq_not_found",
}

_HAPLOTYPE = {"_M": "MATERNAL", "_P": "PATERNAL"}
_CHR_SUFFIX = re.compile(r"(_M|_P)$")


@dataclass
class Normalised:
    rows: list[dict]
    per_read: Path
    summary: Path
    per_arm: Path | None


def _split_chr(value: str) -> tuple[str, str]:
    """`chr01p_M` -> (`chr01p`, `MATERNAL`).

    The R has already truncated MATERNAL to M at posEndPQchrFn; the mapping back is
    unambiguous. Upstream then collapses the two in summary.by_chr; we keep them apart so the
    interface can offer the choice.
    """
    match = _CHR_SUFFIX.search(value)
    if not match:
        return value, ""
    return value[: match.start()], _HAPLOTYPE.get(match.group(1), "")


def read_telomere_r(work: Path) -> list[dict]:
    """Pipeline A. Prefers reads.table.all.txt, which includes the rejected reads."""
    source = work / "r_analysis" / "reads.table.all.txt"
    if not source.is_file():
        source = work / "r_analysis" / "reads.table.txt"
    if not source.is_file():
        raise FileNotFoundError(f"no per-read table under {work / 'r_analysis'}")

    rows = []
    with source.open(newline="") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            # The column is Telomere+Adap upstream, which R's make.names turns into
            # Telomere.Adap when the table is read back. Match on the prefix either way.
            tel_key = next(k for k in record if k.startswith("Telomere"))
            chrom_arm, haplotype = _split_chr(record.get("chr", ""))
            status = record.get("qc_status", "pass")
            length = record.get(tel_key, "")
            rows.append({
                "read_id": record["qname"],
                "telomere_bp": length if status == "pass" else "",
                "read_length_bp": record.get("qwidth", ""),
                "chrom_arm": chrom_arm,
                "haplotype": haplotype,
                "strand": record.get("strand", ""),
                "mapq": record.get("mapq", ""),
                "qc_status": status,
                "qc_reason": record.get("qc_reason", ""),
            })
    return rows


def read_telonp(work: Path) -> list[dict]:
    """Pipeline B. Reference-free, so every mapping column stays empty."""
    candidates = sorted(work.glob("*.telonp.all.csv")) or sorted(work.glob("*.telonp.csv"))
    if not candidates:
        raise FileNotFoundError(f"no *.telonp.all.csv under {work}")

    rows = []
    with candidates[0].open(newline="") as handle:
        for record in csv.DictReader(handle):
            value = int(record["teloNPLength"])
            passed = value > 0
            rows.append({
                "read_id": record["qname"],
                "telomere_bp": str(value) if passed else "",
                "read_length_bp": "",
                "chrom_arm": "",
                "haplotype": "",
                "strand": "",
                "mapq": "",
                "qc_status": "pass" if passed else "fail",
                "qc_reason": "" if passed else TELONP_SENTINELS.get(value, "unknown"),
            })
    return rows


def read_teloseq(work: Path) -> list[dict]:
    """Telo-seq, either variant.

    The pass set is taken from the pipeline's own filtered table rather than re-derived from
    the flags. Re-deriving it is how this first went wrong: the flags alone omit the minimum
    tract length, so the canonical summary disagreed with the pipeline's -- 2,922 "passing"
    reads against the 2,521 the pipeline actually kept. The failure reasons are still worked
    out from the flags, since that is the only place they exist.
    """
    all_csv = sorted(work.glob("*.teloseq.perread.all.csv"))
    if not all_csv:
        raise FileNotFoundError(f"no *.teloseq.perread.all.csv under {work}")
    passed_csv = sorted(work.glob("*.teloseq.perread.csv"))
    arms_csv = sorted(work.glob("*.teloseq.perread.arms.csv"))

    def load(path: Path) -> dict[str, dict]:
        with path.open(newline="") as handle:
            return {r["read_id"]: r for r in csv.DictReader(handle)}

    kept = set(load(passed_csv[0])) if passed_csv else set()
    placed = load(arms_csv[0]) if arms_csv else {}

    rows = []
    with all_csv[0].open(newline="") as handle:
        for record in csv.DictReader(handle):
            read_id = record["read_id"]
            terminal = record.get("is_terminal", "").strip().lower() == "true"
            oriented = record.get("started_correctly", "").strip().lower() == "true"
            arm = (placed.get(read_id) or {}).get("chrom_arm", "")

            sequence_ok = read_id in kept
            passed = sequence_ok and (bool(arm) if arms_csv else True)

            if not terminal:
                reason = "not_terminal"
            elif not oriented:
                reason = "wrong_orientation"
            elif not sequence_ok:
                reason = "tract_too_short"
            elif arms_csv and not arm:
                reason = "no_arm_assignment"
            else:
                reason = ""

            hit = placed.get(read_id) or {}
            rows.append({
                "read_id": read_id,
                "telomere_bp": record.get("telomere_len", "") if passed else "",
                "read_length_bp": record.get("read_len", ""),
                "chrom_arm": arm,
                "haplotype": hit.get("haplotype", ""),
                "strand": "",
                "mapq": hit.get("mapq", ""),
                "qc_status": "pass" if passed else "fail",
                "qc_reason": reason,
            })
    return rows


READERS = {
    "telomere-r": read_telomere_r,
    "telonp": read_telonp,
    "teloseq": read_teloseq,
    "teloseq-mapped": read_teloseq,
}


def _describe(values: list[float]) -> dict:
    if not values:
        return {k: None for k in ("n", "mean", "median", "sd", "min", "max", "q25", "q75")}
    ordered = sorted(values)
    return {
        "n": len(values),
        "mean": round(statistics.fmean(values), 1),
        "median": round(statistics.median(values), 1),
        # R's sd() is the sample standard deviation, so ddof=1 -- statistics.stdev matches.
        "sd": round(statistics.stdev(values), 1) if len(values) > 1 else None,
        "min": int(ordered[0]),
        "max": int(ordered[-1]),
        "q25": round(_quantile(ordered, 0.25), 1),
        "q75": round(_quantile(ordered, 0.75), 1),
    }


def _quantile(ordered: list[float], q: float) -> float:
    """Linear interpolation, matching numpy's default and R's type 7."""
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def summarise(rows: list[dict], *, pipeline_id: str, extra: dict | None = None) -> dict:
    passed = [float(r["telomere_bp"]) for r in rows if r["qc_status"] == "pass" and r["telomere_bp"]]
    breakdown: dict[str, int] = {}
    for row in rows:
        if row["qc_status"] != "pass":
            breakdown[row["qc_reason"] or "unknown"] = breakdown.get(row["qc_reason"] or "unknown", 0) + 1
    return {
        "pipeline": pipeline_id,
        "reads": {
            "total": len(rows),
            "passed": len(passed),
            "pass_rate": round(len(passed) / len(rows), 4) if rows else None,
        },
        "qc_breakdown": dict(sorted(breakdown.items(), key=lambda kv: -kv[1])),
        "telomere_bp": _describe(passed),
        **(extra or {}),
    }


def per_arm(rows: list[dict]) -> list[dict]:
    """Recomputed from the per-read table rather than parsed from the pipeline's own CSV, so a
    filter applied later (a mapq cutoff, say) gives correct numbers without re-running."""
    groups: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        if row["qc_status"] != "pass" or not row["chrom_arm"] or not row["telomere_bp"]:
            continue
        groups.setdefault((row["chrom_arm"], row["haplotype"]), []).append(float(row["telomere_bp"]))
    out = []
    for (arm, haplotype), values in sorted(groups.items()):
        stats = _describe(values)
        out.append({
            "chrom_arm": arm, "haplotype": haplotype, "n": stats["n"],
            "median_telomere_bp": stats["median"], "mean_telomere_bp": stats["mean"],
            "sd_telomere_bp": stats["sd"],
        })
    return out


def normalise(work: Path, out: Path, *, pipeline_id: str, extra: dict | None = None) -> Normalised:
    """Read a finished run's outputs and write the canonical files."""
    rows = READERS[pipeline_id](work)
    out.mkdir(parents=True, exist_ok=True)

    per_read_path = out / "per_read.csv"
    with per_read_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PER_READ_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    arms = per_arm(rows)
    per_arm_path = None
    if arms:
        per_arm_path = out / "per_arm.csv"
        with per_arm_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(arms[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(arms)

    summary_path = out / "summary.json"
    summary_path.write_text(
        json.dumps(summarise(rows, pipeline_id=pipeline_id, extra=extra), indent=2) + "\n"
    )
    return Normalised(rows, per_read_path, summary_path, per_arm_path)
