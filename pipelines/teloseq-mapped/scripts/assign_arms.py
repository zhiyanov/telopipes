#!/usr/bin/env python3
"""Attach a chromosome arm to each Telo-seq read.

NCRF measures the telomere from the read's own sequence and has no idea which chromosome the
read came from. Mapping supplies that, and only that: the length still comes from NCRF.

**This is not Telo-seq's own anchoring.** The published pipeline maps to a reference that
retains its telomeres and checks, in filter_bam_2023.py, that the alignment spans telomere and
subtelomere intervals listed in a per-assembly anchor TSV. We instead map to the same cut
reference the mapping-based pipeline uses -- telomeres removed, one record per arm -- and
require the alignment to start close to the subtelomere boundary. The intent is the same (the
read must genuinely reach a chromosome end) but the mechanism differs, so arm assignments need
not agree read-for-read with the published tool. Recorded in docs/PATCHES.md.

Using the cut reference is what keeps a single reference type in the application: the same
registered reference serves this pipeline and the mapping-based one.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys

CHR_ARM = re.compile(r"^(chr(?:[0-9]{2}|X|Y)[pq])(?:_(MATERNAL|PATERNAL))?$")

FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


def primary_alignments(bam: str, arm_length: int, max_map_start: int):
    """Yield read_id -> (arm, haplotype, mapq, distance from the subtelomere boundary).

    A p-arm record starts at the boundary, so distance is the leftmost reference position. A
    q-arm record *ends* at it, so distance is measured back from the record length -- the same
    mirroring the R pipeline does by flipping q-arm coordinates negative.
    """
    proc = subprocess.Popen(
        ["samtools", "view", "-F", str(FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_UNMAPPED), bam],
        stdout=subprocess.PIPE, text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.split("\t", 9)
        if len(fields) < 9:
            continue
        qname, _, rname, pos, mapq, cigar = fields[0], fields[1], fields[2], int(fields[3]), int(fields[4]), fields[5]
        match = CHR_ARM.match(rname)
        if not match:
            continue
        arm, haplotype = match.group(1), match.group(2) or ""

        if arm.endswith("p"):
            distance = pos
        else:
            span = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHPX=])", cigar) if op in "MDN=X")
            distance = arm_length - (pos + span - 1)
        if distance <= max_map_start:
            yield qname, (arm, haplotype, mapq, max(distance, 0))
    proc.wait()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("perread", help="<sample>.teloseq.perread.csv from the length step")
    ap.add_argument("bam", help="coordinate-sorted BAM against the cut reference")
    ap.add_argument("out", help="output CSV with arm assignment")
    ap.add_argument("--arm-length", type=int, required=True)
    ap.add_argument("--map-start-threshold", type=int, default=1000,
                    help="max distance of the alignment start from the subtelomere boundary "
                         "[%(default)s]")
    args = ap.parse_args()

    placed = dict(primary_alignments(args.bam, args.arm_length, args.map_start_threshold))

    with open(args.perread, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        print("error: no reads in the per-read table", file=sys.stderr)
        return 1

    fields = list(rows[0]) + ["chrom_arm", "haplotype", "mapq", "map_start_from_boundary"]
    assigned = 0
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            hit = placed.get(row["read_id"])
            row.update(zip(fields[-4:], hit if hit else ("", "", "", "")))
            assigned += hit is not None
            writer.writerow(row)

    print(f"[teloseq] {assigned:,} of {len(rows):,} reads assigned to a chromosome arm "
          f"({assigned / len(rows):.1%})")
    print(f"[teloseq] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
