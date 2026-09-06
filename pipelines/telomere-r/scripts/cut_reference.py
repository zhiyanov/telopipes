#!/usr/bin/env python3
"""Derive a cut reference from a whole genome, for the mapping-based pipeline.

`telomere-r` measures a telomere as the bases minimap2 had to soft-clip. That only works if
the reference has no telomere to align against: give it a telomere-to-telomere assembly and
the aligner extends into the assembly's own telomere, so the clip becomes
(sample telomere - assembly telomere) and any read shorter than the assembly's telomere
registers as nothing at all. So the telomere must be removed, and each chromosome end kept as
its own fixed-length record.

This script does that, and normalises the naming while it is at it:

    chr1_MATERNAL (244 Mbp)  ->  chr01p_MATERNAL (500 kb)   from the p end, telomere trimmed
                             ->  chr01q_MATERNAL (500 kb)   from the q end, telomere trimmed

Normalising here rather than trusting the downstream R matters. posEndPQchrFn() zero-pads
chr1p to chr01p on one line and then overwrites the result from rname on the next, so its own
repair is dead code -- and an unpadded name does not error, it just sorts wrongly everywhere.
Emitting canonical names is the only place that gets fixed.

Deliberately does not use TeloBP's trimTeloReferenceGenome(): that emits _C/_G suffixes rather
than p/q, and when a record's end telomere length is zero it evaluates record[-0-N:-0], i.e.
record[-N:0], which is empty.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

#: Telomeric hexamer and its common variant repeats, on both strands. p ends carry the C-rich
#: strand, q ends the G-rich one. Deliberately no five-mers: CCTAA alone occurs by chance
#: every ~1 kb of random sequence, which is enough to fake a short tract at a telomere-less end.
P_REPEAT = re.compile(r"CCCTAA|CCCCTAA|CCCTCA|CCCTGA")
Q_REPEAT = re.compile(r"TTAGGG|TTAGGGG|TGAGGG|TCAGGG")

#: Gap tolerated inside a telomere tract before it is treated as ended. 250 bp mirrors
#: Telo-seq's build_telomere_reference.py, which merges repeat intervals with the same limit.
MAX_GAP = 250

#: A candidate tract must be at least this long and this densely matched to count. Without a
#: floor, a couple of chance hexamers near a telomere-less terminus trim a spurious ~150 bp --
#: which silently shifts that arm's coordinate frame relative to every other arm.
MIN_TRACT_BP = 100
MIN_TRACT_DENSITY = 0.6

#: Records we never cut: no telomeres, or not a chromosome.
SKIP = re.compile(
    r"(^|_)(M|MT|chrM|chrMT|EBV)($|_)|mitochond|scaffold|unplaced|unlocalized|random|_alt|decoy",
    re.IGNORECASE,
)

CHROM = re.compile(r"(?:^|[^0-9A-Za-z])(?:chr)?([0-9]{1,2}|X|Y)(?:$|[^0-9A-Za-z])", re.IGNORECASE)
HAPLOTYPE = [
    (re.compile(r"MATERNAL|\bmat\b|hap1|haplotype1|\.1$", re.IGNORECASE), "MATERNAL"),
    (re.compile(r"PATERNAL|\bpat\b|hap2|haplotype2|\.2$", re.IGNORECASE), "PATERNAL"),
]


@dataclass
class ArmReport:
    name: str
    source: str
    arm: str
    telomere_bp: int
    source_length: int


@dataclass
class Report:
    arm_length_bp: int
    arms: list[ArmReport] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "arm_length_bp": self.arm_length_bp,
            "n_arms": len(self.arms),
            "haplotypes": sorted({a.name.split("_")[1] for a in self.arms if "_" in a.name}),
            "arms": [vars(a) for a in self.arms],
            "skipped": [{"record": r, "reason": why} for r, why in self.skipped],
        }


def canonical_name(header: str) -> tuple[str, str] | tuple[None, str]:
    """`chr1_MATERNAL` -> (`chr01`, `MATERNAL`). Returns (None, reason) if unusable."""
    record = header.split()[0]
    if SKIP.search(record):
        return None, "not a chromosome, or has no telomeres"
    match = CHROM.search(record)
    if not match:
        return None, "no chromosome number found in the name"
    token = match.group(1).upper()
    chrom = f"chr{int(token):02d}" if token.isdigit() else f"chr{token}"
    for pattern, label in HAPLOTYPE:
        if pattern.search(record):
            return chrom, label
    return chrom, ""


def telomere_length(seq: str, *, at_start: bool, search: int) -> int:
    """Length of the telomeric tract at one end of a record.

    Repeat matches are walked outward from the terminus and the tract ends at the first gap
    wider than MAX_GAP. Returns 0 when the end carries no telomere -- which is a real case for
    an incomplete assembly, and the caller must handle it rather than slicing by a zero offset.
    """
    window = seq[:search] if at_start else seq[-search:]
    pattern = P_REPEAT if at_start else Q_REPEAT
    spans = [(m.start(), m.end()) for m in pattern.finditer(window)]
    if not spans:
        return 0

    if at_start:
        if spans[0][0] > MAX_GAP:
            return 0                      # nothing telomeric near the terminus
        end, used = spans[0][1], [spans[0]]
        for start, stop in spans[1:]:
            if start - end > MAX_GAP:
                break
            end, _ = stop, used.append((start, stop))
        tract = end
    else:
        if len(window) - spans[-1][1] > MAX_GAP:
            return 0
        start, used = spans[-1][0], [spans[-1]]
        for span_start, span_stop in reversed(spans[:-1]):
            if start - span_stop > MAX_GAP:
                break
            start, _ = span_start, used.append((span_start, span_stop))
        tract = len(window) - start

    matched = sum(stop - begin for begin, stop in used)
    if tract < MIN_TRACT_BP or matched / tract < MIN_TRACT_DENSITY:
        return 0                          # chance matches, not a telomere
    return tract


def iter_records(path: Path, keep: int):
    """Stream (header, head, tail, total_length).

    Only the first and last `keep` bases of each record are retained, so a 6 Gbp diploid
    assembly costs a couple of megabytes rather than being read into memory whole.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    header, head, tail, total = None, [], [], 0
    head_len = 0

    def flush():
        return header, "".join(head), "".join(tail)[-keep:], total

    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield flush()
                header, head, tail, total, head_len = line[1:].strip(), [], [], 0, 0
                continue
            if header is None:
                continue
            seq = line.strip()
            total += len(seq)
            if head_len < keep:
                head.append(seq)
                head_len += len(seq)
            tail.append(seq)
            if len(tail) > 64:                     # bound the tail buffer
                joined = "".join(tail)
                tail = [joined[-keep:]] if len(joined) > keep else [joined]
    if header is not None:
        yield flush()


def write_record(handle, name: str, seq: str, width: int = 60) -> None:
    handle.write(f">{name}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def cut(source: Path, out_fasta: Path, *, arm_length: int, search: int,
        haplotype: str = "all") -> Report:
    report = Report(arm_length_bp=arm_length)
    keep = arm_length + search

    with out_fasta.open("w") as out:
        for header, head, tail, total in iter_records(source, keep):
            chrom, hap_or_reason = canonical_name(header)
            record = header.split()[0]
            if chrom is None:
                report.skipped.append((record, hap_or_reason))
                continue
            hap = hap_or_reason
            if haplotype != "all" and hap and hap != haplotype:
                report.skipped.append((record, f"not the requested {haplotype} haplotype"))
                continue
            suffix = f"_{hap}" if hap else ""

            if total < 2 * arm_length:
                report.skipped.append(
                    (record, f"only {total:,} bp; needs at least {2 * arm_length:,} "
                             f"for two non-overlapping {arm_length:,} bp arms")
                )
                continue

            for arm, at_start, chunk in (("p", True, head), ("q", False, tail)):
                telo = telomere_length(chunk, at_start=at_start, search=search)
                # p arm runs inward from just past the telomere; q arm runs inward from just
                # before it. Slicing explicitly, because the zero-telomere case is exactly
                # where upstream's helper produces an empty record.
                body = chunk[telo : telo + arm_length] if at_start else (
                    chunk[len(chunk) - telo - arm_length : len(chunk) - telo] if telo
                    else chunk[-arm_length:]
                )
                if len(body) != arm_length:
                    report.skipped.append(
                        (f"{record} {arm} arm",
                         f"only {len(body):,} bp available after trimming {telo:,} bp of telomere")
                    )
                    continue
                name = f"{chrom}{arm}{suffix}"
                write_record(out, name, body)
                report.arms.append(ArmReport(name, record, arm, telo, total))
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome", type=Path, help="whole-genome FASTA, optionally gzipped")
    ap.add_argument("out", type=Path, help="output cut reference FASTA")
    ap.add_argument("--arm-length", type=int, default=500_000,
                    help="bases kept per chromosome end [%(default)s]")
    ap.add_argument("--telomere-search", type=int, default=200_000,
                    help="bases searched at each end for the telomere tract [%(default)s]")
    ap.add_argument("--haplotype", choices=["MATERNAL", "PATERNAL", "all"], default="all",
                    help="keep only one haplotype of a diploid assembly [%(default)s]")
    ap.add_argument("--report", type=Path, help="write a JSON report of what was cut")
    args = ap.parse_args()

    if not args.genome.is_file():
        print(f"error: no such file: {args.genome}", file=sys.stderr)
        return 2

    args.out.parent.mkdir(parents=True, exist_ok=True)
    report = cut(args.genome, args.out, arm_length=args.arm_length,
                 search=args.telomere_search, haplotype=args.haplotype)

    if not report.arms:
        print("error: no usable chromosome arms were produced", file=sys.stderr)
        for record, why in report.skipped[:10]:
            print(f"  skipped {record}: {why}", file=sys.stderr)
        return 1

    if args.report:
        args.report.write_text(json.dumps(report.as_dict(), indent=2) + "\n")

    haplotypes = sorted({a.name.split("_")[1] for a in report.arms if "_" in a.name})
    telos = [a.telomere_bp for a in report.arms]
    print(f"[cut] {len(report.arms)} arms x {args.arm_length:,} bp -> {args.out}")
    print(f"[cut] haplotypes: {', '.join(haplotypes) or 'none (haploid, unlabelled)'}")
    print(f"[cut] telomere trimmed per arm: median {sorted(telos)[len(telos) // 2]:,} bp, "
          f"range {min(telos):,}-{max(telos):,} bp")
    if any(t == 0 for t in telos):
        n = sum(1 for t in telos if t == 0)
        print(f"[cut] note: {n} arm(s) had no detectable telomere -- the assembly does not "
              f"reach the end there, so reads will clip against an internal boundary")
    noise = sum(1 for _, why in report.skipped if "haplotype" in why)
    for record, why in report.skipped:
        if "haplotype" not in why:
            print(f"[cut] skipped {record}: {why}")
    if noise:
        print(f"[cut] skipped {noise} records from the other haplotype")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
