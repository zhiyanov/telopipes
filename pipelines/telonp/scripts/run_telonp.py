#!/usr/bin/env python3
"""Method 2 (sequence-based): per-read telomere length via GreiderLab TeloNP.

Runs getTeloNPBoundary() on every read of an ONT FASTQ and writes a CSV of per-read telomere
lengths plus a summary and histogram. TeloNP is TeloBP tuned for ONT Guppy v6 super-high-accuracy
basecaller error patterns (Karimian et al., Science 2024). Reads that fail QC return negative
sentinels and are dropped, the same convention as the upstream teloBPCmd.py.

Differences from the reference driver, all recorded in docs/PATCHES.md:

  * Reads are streamed and scored in bounded chunks via ProcessPoolExecutor rather than loaded
    whole into a DataFrame and dispatched with pandarallel. The reference's peak memory scaled
    as (total sequence bytes x workers) because pandarallel pickles a DataFrame slice per
    worker through /dev/shm -- which podman caps at 64 MB, so the reference driver fails
    outright in a container with ENOSPC unless --shm-size is raised. Peak memory here is
    O(chunk_reads x workers) and no shared memory is involved.
  * pandas is gone from the compute path entirely, so the pandas version no longer affects
    anything.
  * .fastq.gz is accepted; SeqIO.parse(path, "fastq") does not decompress.
  * The sys.path surgery the reference needed is gone. It existed because the vendored TeloBP
    checkout sat next to the driver and shadowed the installed package as a namespace package;
    in the image TeloBP is in site-packages and this file is in /opt/telomers, so that cannot
    happen. A failed import is reported explicitly instead of silently half-working.

Output is byte-identical to the reference driver's.
"""
import argparse
import csv
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import islice

from Bio import SeqIO

try:
    from TeloBP import getTeloNPBoundary
except ImportError as exc:  # pragma: no cover - environment problem, not logic
    sys.exit(
        f"could not import getTeloNPBoundary from TeloBP: {exc}\n"
        "If a TeloBP source checkout sits next to this script it can shadow the installed "
        "package as a namespace package. Run from a directory that does not contain one."
    )

# Same sentinel map the upstream package uses for failed reads.
ERROR_RETURNS = {"init": -1, "fusedRead": -10, "strandType": -20, "seqNotFound": -1000}


def score(seq):
    """Worker entry point. Must be module-level to be picklable."""
    if not seq:
        return ERROR_RETURNS["seqNotFound"]
    return getTeloNPBoundary(seq)


def iter_reads(path):
    """Stream (id, sequence) pairs, transparently handling gzip."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            yield record.id, str(record.seq)


def chunked(iterable, size):
    it = iter(iterable)
    while chunk := list(islice(it, size)):
        yield chunk


def write_histogram(lengths, sample, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    arr = np.asarray(lengths)
    plt.figure(figsize=(8, 5))
    upper = np.quantile(arr, 0.99)
    plt.hist(arr, bins=np.arange(0, upper + 250, 250), color="#2c7fb8", edgecolor="white")
    plt.axvline(np.median(arr), color="crimson", lw=2, label=f"median {np.median(arr):.0f} bp")
    plt.xlabel("Telomere length (bp)")
    plt.ylabel("Reads")
    plt.title(f"TeloNP telomere length — {sample} (n={len(arr)})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=120)
    plt.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fastq", help="input FASTQ (optionally gzipped)")
    ap.add_argument("outdir", help="output directory")
    ap.add_argument("--sample", default=None, help="sample label (default: fastq basename)")
    ap.add_argument("--workers", type=int, default=0, help="parallel workers (0 = all cores)")
    ap.add_argument("--chunk-reads", type=int, default=2000,
                    help="reads held in memory per dispatch round; bounds peak RSS")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    name = os.path.basename(args.fastq)
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    sample = args.sample or name.rsplit(".", 1)[0]
    workers = args.workers if args.workers > 0 else (os.cpu_count() or 1)

    all_csv = os.path.join(args.outdir, f"{sample}.telonp.all.csv")
    passed_csv = os.path.join(args.outdir, f"{sample}.telonp.csv")

    print(f"[telonp] streaming reads from {args.fastq} ({workers} workers)")

    total = 0
    lengths = []                                    # positive lengths only
    counts = {name: 0 for name in ERROR_RETURNS}
    by_value = {value: name for name, value in ERROR_RETURNS.items()}

    with open(all_csv, "w", newline="") as fh_all, \
         open(passed_csv, "w", newline="") as fh_pass, \
         ProcessPoolExecutor(max_workers=workers) as pool:
        w_all = csv.writer(fh_all, lineterminator="\n")
        w_pass = csv.writer(fh_pass, lineterminator="\n")
        w_all.writerow(["qname", "teloNPLength"])
        w_pass.writerow(["qname", "teloNPLength"])

        for chunk in chunked(iter_reads(args.fastq), args.chunk_reads):
            qnames = [qname for qname, _ in chunk]
            seqs = [seq for _, seq in chunk]
            for qname, length in zip(qnames, pool.map(score, seqs)):
                total += 1
                w_all.writerow([qname, length])
                if length > 0:
                    lengths.append(length)
                    w_pass.writerow([qname, length])
                elif length in by_value:
                    counts[by_value[length]] += 1
            print(f"[telonp] {total} reads scored", flush=True)

    n_pass = len(lengths)
    if n_pass:
        import numpy as np
        arr = np.asarray(lengths)
        stats = [
            f"mean_telomere_bp\t{arr.mean():.1f}",
            f"median_telomere_bp\t{np.median(arr):.1f}",
            f"sd_telomere_bp\t{arr.std(ddof=1):.1f}" if n_pass > 1 else "sd_telomere_bp\tNA",
            f"min_telomere_bp\t{int(arr.min())}",
            f"max_telomere_bp\t{int(arr.max())}",
        ]
    else:
        stats = [f"{k}_telomere_bp\tNA" for k in ("mean", "median", "sd", "min", "max")]

    summary_lines = [
        f"sample\t{sample}",
        f"total_reads\t{total}",
        f"passed_reads\t{n_pass}",
        f"pass_rate\t{n_pass / max(total, 1):.4f}",
        *stats,
        *(f"err_{name}\t{counts[name]}" for name in ERROR_RETURNS),
    ]
    summary = "\n".join(summary_lines) + "\n"
    with open(os.path.join(args.outdir, f"{sample}.telonp.summary.txt"), "w") as fh:
        fh.write(summary)
    print("[telonp] summary:\n" + summary)

    if n_pass:
        write_histogram(lengths, sample, os.path.join(args.outdir, f"{sample}.telonp.hist.png"))

    print(f"[telonp] wrote {passed_csv}")


if __name__ == "__main__":
    main()
