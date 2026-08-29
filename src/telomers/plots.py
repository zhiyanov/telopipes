"""Figures rendered from the canonical tables.

These sit alongside, not instead of, the figures each pipeline emits. The pipelines' own PNGs
are the citable ones and are kept as artifacts -- but pipeline A produces ~19 of them in its
own style with four hardcoded mapq cutoffs, while pipeline B produces one, so a consistent
interface needs figures drawn from the common schema.
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK = "#2c7fb8"
ACCENT = "#c1352a"


def _passed(rows: list[dict]) -> list[float]:
    return [float(r["telomere_bp"]) for r in rows if r["qc_status"] == "pass" and r["telomere_bp"]]


def read_per_read(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def distribution(rows: list[dict], out: Path, *, sample: str = "") -> Path | None:
    values = _passed(rows)
    if not values:
        return None
    ordered = sorted(values)
    upper = ordered[int(len(ordered) * 0.99)]
    median = statistics.median(values)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = [i for i in range(0, int(upper) + 250, 250)]
    ax.hist(values, bins=bins, color=INK, edgecolor="white", linewidth=0.5)
    ax.axvline(median, color=ACCENT, lw=2, label=f"median {median:,.0f} bp")
    ax.set_xlabel("Telomere length (bp)")
    ax.set_ylabel("Reads")
    ax.set_title(f"Telomere length distribution — {sample} (n={len(values):,})".strip(" —"))
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def by_arm(rows: list[dict], out: Path, *, sample: str = "") -> Path | None:
    """Per chromosome arm, ordered chr01p ... chrXq. Only meaningful for pipeline A."""
    groups: dict[str, list[float]] = {}
    for row in rows:
        if row["qc_status"] != "pass" or not row["chrom_arm"] or not row["telomere_bp"]:
            continue
        groups.setdefault(row["chrom_arm"], []).append(float(row["telomere_bp"]))
    if not groups:
        return None

    arms = sorted(groups)
    fig, ax = plt.subplots(figsize=(max(8, len(arms) * 0.28), 5))
    ax.boxplot([groups[a] for a in arms], tick_labels=arms, showfliers=False,
               medianprops={"color": ACCENT, "linewidth": 1.5},
               boxprops={"color": "#555"}, whiskerprops={"color": "#555"},
               capprops={"color": "#555"})
    overall = statistics.median([v for values in groups.values() for v in values])
    ax.axhline(overall, color=INK, lw=1, ls="--", alpha=0.7,
               label=f"overall median {overall:,.0f} bp")
    ax.set_ylabel("Telomere length (bp)")
    ax.set_title(f"Telomere length by chromosome arm — {sample}".strip(" —"))
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def read_funnel(rows: list[dict], out: Path, *, sample: str = "") -> Path | None:
    """How many reads survived, and why the rest did not.

    Worth drawing: pipeline A discarded this information entirely until we added --emit-all,
    so a low yield used to be invisible.
    """
    if not rows:
        return None
    reasons: dict[str, int] = {}
    for row in rows:
        if row["qc_status"] != "pass":
            reasons[row["qc_reason"] or "unknown"] = reasons.get(row["qc_reason"] or "unknown", 0) + 1
    passed = sum(1 for r in rows if r["qc_status"] == "pass")

    labels = ["passed"] + sorted(reasons, key=lambda k: -reasons[k])
    counts = [passed] + [reasons[k] for k in labels[1:]]
    colours = [INK] + ["#b9b9b9"] * (len(labels) - 1)

    fig, ax = plt.subplots(figsize=(7, max(2.2, 0.5 * len(labels) + 1)))
    positions = range(len(labels))
    ax.barh(list(positions), counts, color=colours)
    ax.set_yticks(list(positions), [k.replace("_", " ") for k in labels])
    ax.invert_yaxis()
    for y, count in zip(positions, counts):
        ax.text(count, y, f" {count:,}", va="center", fontsize=9)
    ax.set_xlabel("Reads")
    ax.set_title(f"Read accounting — {sample} ({len(rows):,} total)".strip(" —"))
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, max(counts) * 1.15)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def render_all(per_read: Path, out_dir: Path, *, sample: str = "") -> list[Path]:
    rows = read_per_read(per_read)
    out_dir.mkdir(parents=True, exist_ok=True)
    produced = [
        distribution(rows, out_dir / "distribution.png", sample=sample),
        by_arm(rows, out_dir / "by_arm.png", sample=sample),
        read_funnel(rows, out_dir / "read_funnel.png", sample=sample),
    ]
    return [p for p in produced if p is not None]
