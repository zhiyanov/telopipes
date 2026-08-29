#!/usr/bin/env python3
"""Method 2 (sequence-based): per-read telomere length via GreiderLab TeloNP.

Runs getTeloNPBoundary() on every read of an ONT FASTQ and writes a CSV of
per-read telomere lengths plus a summary and histogram. TeloNP is TeloBP tuned
for ONT Guppy v6 super-high-accuracy basecaller error patterns (Karimian et al.,
Science 2024). Reads that fail QC return negative sentinels and are dropped, the
same convention as the upstream teloBPCmd.py.
"""
import argparse
import os
import sys

# The vendored TeloBP repo lives next to this script (pipeline/TeloBP/), whose
# outer folder would shadow the pip-installed package as a namespace package.
# Drop this script's own directory from sys.path so the real install wins.
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:] = [p for p in sys.path if os.path.abspath(p or ".") != _here]

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from pandarallel import pandarallel

from TeloBP import getTeloNPBoundary

# Same sentinel map the upstream package uses for failed reads.
ERROR_RETURNS = {"init": -1, "fusedRead": -10, "strandType": -20, "seqNotFound": -1000}


def row_to_telonp(row):
    seq = row["seq"]
    if seq is np.nan:
        return ERROR_RETURNS["seqNotFound"]
    return getTeloNPBoundary(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fastq", help="input FASTQ")
    ap.add_argument("outdir", help="output directory")
    ap.add_argument("--sample", default=None, help="sample label (default: fastq basename)")
    ap.add_argument("--workers", type=int, default=0, help="parallel workers (0 = all cores)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sample = args.sample or os.path.basename(args.fastq).rsplit(".", 1)[0]

    print(f"[telonp] loading reads from {args.fastq}")
    records = [(r.id, str(r.seq)) for r in SeqIO.parse(args.fastq, "fastq")]
    df = pd.DataFrame(records, columns=["qname", "seq"])
    print(f"[telonp] {len(df)} reads loaded; scoring with TeloNP...")

    if args.workers and args.workers > 0:
        pandarallel.initialize(nb_workers=args.workers, progress_bar=True)
    else:
        pandarallel.initialize(progress_bar=True)
    df["teloNPLength"] = df.parallel_apply(row_to_telonp, axis=1)

    df = df.drop(columns=["seq"])
    all_csv = os.path.join(args.outdir, f"{sample}.telonp.all.csv")
    df.to_csv(all_csv, index=False)

    passed = df[df["teloNPLength"] > 0].reset_index(drop=True)
    passed_csv = os.path.join(args.outdir, f"{sample}.telonp.csv")
    passed.to_csv(passed_csv, index=False)

    # Error breakdown by sentinel.
    counts = {name: int((df["teloNPLength"] == val).sum()) for name, val in ERROR_RETURNS.items()}
    lengths = passed["teloNPLength"].to_numpy()

    summary_lines = [
        f"sample\t{sample}",
        f"total_reads\t{len(df)}",
        f"passed_reads\t{len(passed)}",
        f"pass_rate\t{len(passed) / max(len(df), 1):.4f}",
        f"mean_telomere_bp\t{lengths.mean():.1f}" if len(lengths) else "mean_telomere_bp\tNA",
        f"median_telomere_bp\t{np.median(lengths):.1f}" if len(lengths) else "median_telomere_bp\tNA",
        f"sd_telomere_bp\t{lengths.std(ddof=1):.1f}" if len(lengths) > 1 else "sd_telomere_bp\tNA",
        f"min_telomere_bp\t{int(lengths.min())}" if len(lengths) else "min_telomere_bp\tNA",
        f"max_telomere_bp\t{int(lengths.max())}" if len(lengths) else "max_telomere_bp\tNA",
    ]
    for name, val in ERROR_RETURNS.items():
        summary_lines.append(f"err_{name}\t{counts[name]}")
    summary = "\n".join(summary_lines) + "\n"
    with open(os.path.join(args.outdir, f"{sample}.telonp.summary.txt"), "w") as fh:
        fh.write(summary)
    print("[telonp] summary:\n" + summary)

    if len(lengths):
        plt.figure(figsize=(8, 5))
        upper = np.quantile(lengths, 0.99)
        plt.hist(lengths, bins=np.arange(0, upper + 250, 250), color="#2c7fb8", edgecolor="white")
        plt.axvline(np.median(lengths), color="crimson", lw=2,
                    label=f"median {np.median(lengths):.0f} bp")
        plt.xlabel("Telomere length (bp)")
        plt.ylabel("Reads")
        plt.title(f"TeloNP telomere length — {sample} (n={len(lengths)})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"{sample}.telonp.hist.png"), dpi=120)
        plt.close()

    print(f"[telonp] wrote {passed_csv}")


if __name__ == "__main__":
    main()
