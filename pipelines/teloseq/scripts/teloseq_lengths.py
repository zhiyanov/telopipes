#!/usr/bin/env python3
"""Telo-seq step 2: per-read telomere length from NCRF's alignment summary.

Ports the three functions the reference implementation loaded out of Telo-seq's
`scripts/ncrf.summary.stats.py` with importlib. Reimplemented here so the pipeline does not
depend on a path inside a vendored 20-conda-environment Snakemake repo of which one file was
ever used, and so the whole thing runs inside the image.

What the length actually is, and why it differs from the other pipelines:

    telomere_len = end.max() - start.min()

i.e. the span from the first to the last telomeric repeat NCRF found in the read. If a read
carries an interstitial telomeric repeat well inside the subtelomere, that span reaches across
the intervening non-telomeric sequence and the length is inflated. `segments` records how many
separate NCRF alignments contributed, which is the diagnostic for exactly that case. TeloNP and
the mapping pipeline both measure to a *boundary* instead, so their numbers are not the same
quantity even though the bulk medians land close together.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

COMPLEMENT = str.maketrans("ATGC", "TACG")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def aggregate(group: pd.DataFrame) -> pd.Series:
    """One row per read, from however many NCRF alignments it produced."""
    start, end, seq_len = group["start"].min(), group["end"].max(), group["seqLen"].max()
    motif = group["motif"].max()
    values = {
        "mstart": start,
        "mend": end,
        "telomere_len": end - start,
        "motif": revcomp(motif) if group["strand"].max() == "-" else motif,
        "subtelomere_len": seq_len - (end - start),
        "read_len": seq_len,
        "segments": len(group),
    }
    for column in ("m", "mm", "i", "d"):
        values[column] = group[column].sum()
    return pd.Series(values)


def started_correctly(row, padding: int = 200) -> bool:
    """A genuine terminal telomere reads CCCTAA at the read start (p arm) or TTAGGG at the
    read end (q arm). Anything else is repeat sitting somewhere it should not be."""
    if row.motif == "CCCTAA":
        return bool(row.mstart < padding)
    if row.motif == "TTAGGG":
        return bool(row.mend > row.read_len - padding)
    return False


def is_terminal(row, padding: int = 200) -> bool:
    return bool(row.mstart < padding or row.mend > row.read_len - padding)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("summary", help="reads.ncrf.summary from NCRF")
    ap.add_argument("outdir")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-telomere", type=int, default=1000,
                    help="drop reads whose tract is shorter than this [%(default)s]")
    ap.add_argument("--padding", type=int, default=200,
                    help="how close to a read end a tract must start/finish [%(default)s]")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.summary, sep="\t",
                     usecols=["seq", "start", "end", "seqLen", "motif", "strand",
                              "m", "mm", "i", "d"])
    df = df.rename(columns={"seq": "read_id"})
    df["motif"] = df["motif"].str.replace("telomeric", "TTAGGG")

    per = df.groupby("read_id").apply(aggregate, include_groups=False)
    per["ident"] = 100 * (per.m / per[["m", "mm", "i", "d"]].sum(axis=1))
    per["has_subtelomere"] = per.subtelomere_len > 500
    per["started_correctly"] = per.apply(started_correctly, axis=1, padding=args.padding)
    per["is_terminal"] = per.apply(is_terminal, axis=1, padding=args.padding)
    per = per.reset_index()

    all_csv = os.path.join(args.outdir, f"{args.sample}.teloseq.perread.all.csv")
    per.to_csv(all_csv, index=False)

    keep = per[per.is_terminal & per.started_correctly.fillna(False)
               & (per.telomere_len >= args.min_telomere)].reset_index(drop=True)
    passed_csv = os.path.join(args.outdir, f"{args.sample}.teloseq.perread.csv")
    keep.to_csv(passed_csv, index=False)

    lengths = keep.telomere_len.to_numpy(dtype=float)
    lines = [
        f"sample\t{args.sample}",
        f"reads_with_telomere_tract\t{len(per)}",
        f"reads_passing_filters\t{len(keep)}",
        f"min_telomere_filter_bp\t{args.min_telomere}",
    ]
    if len(lengths):
        lines += [
            f"mean_telomere_bp\t{lengths.mean():.1f}",
            f"median_telomere_bp\t{np.median(lengths):.1f}",
            f"sd_telomere_bp\t{lengths.std(ddof=1):.1f}" if len(lengths) > 1 else "sd_telomere_bp\tNA",
            f"min_telomere_bp\t{int(lengths.min())}",
            f"max_telomere_bp\t{int(lengths.max())}",
        ]
    else:
        lines += [f"{k}_telomere_bp\tNA" for k in ("mean", "median", "sd", "min", "max")]
    for motif, group in keep.groupby("motif"):
        lines.append(f"n_{motif}\t{len(group)}\tmedian_{motif}_bp\t{group.telomere_len.median():.0f}")

    summary = "\n".join(lines) + "\n"
    with open(os.path.join(args.outdir, f"{args.sample}.teloseq.summary.txt"), "w") as fh:
        fh.write(summary)
    print("[teloseq] summary:\n" + summary)

    if len(lengths):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 5))
        upper = np.quantile(lengths, 0.99)
        plt.hist(lengths, bins=np.arange(0, upper + 250, 250), color="#2c7fb8", edgecolor="white")
        plt.axvline(np.median(lengths), color="crimson", lw=2,
                    label=f"median {np.median(lengths):.0f} bp")
        plt.xlabel("Telomere tract length (bp)")
        plt.ylabel("Reads")
        plt.title(f"Telo-seq / NCRF tract length — {args.sample} (n={len(lengths)})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"{args.sample}.teloseq.hist.png"), dpi=120)
        plt.close()

    print(f"[teloseq] wrote {passed_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
