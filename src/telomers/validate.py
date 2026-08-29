"""Input validation.

The reference validator is the most important code in this module, and possibly in the project.
Pipeline A's q-arm telomere length is computed as `softClip - |posEnd| + 1` where `posEnd` was
derived by subtracting `chrArmLn`. If that value disagrees with the reference's actual record
length, **q-arm lengths are silently corrupted while p arms stay perfect** -- a half-broken
result that looks like biology rather than a bug. So the arm length is derived from the file
and never accepted from the user.

Contig naming is load-bearing too, and upstream's own repair for it is dead code:
`posEndPQchrFn` zero-pads `chr1p` to `chr01p` on one line and then overwrites `reads$chr` from
`reads$rname` on the next, discarding it. Unpadded names do not error -- the arm regex still
matches -- they just produce a differently sorted label set. Catching that here is the only
place it gets caught.
"""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

#: chr01p_MATERNAL, chr22q, chrXp_PATERNAL. Zero-padded autosomes only, because the R will not
#: pad them for us and unpadded names sort wrongly in every per-chromosome output.
CONTIG_RE = re.compile(r"^chr(0[1-9]|1[0-9]|2[0-2]|X|Y)([pq])(_(MATERNAL|PATERNAL))?$")


@dataclass
class ReferenceInfo:
    n_records: int = 0
    arm_length_bp: int | None = None
    contigs: list[str] = field(default_factory=list)
    haplotypes: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def _iter_fasta_lengths(path: Path):
    """Yield (header, length) without holding any sequence in memory."""
    name, length = None, 0
    with _open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, length
                name, length = line[1:].strip().split()[0] if line[1:].strip() else "", 0
            else:
                length += len(line.strip())
    if name is not None:
        yield name, length


def validate_reference(path: Path, *, max_report: int = 8) -> ReferenceInfo:
    """Check a cut reference and derive its arm length.

    A "cut" reference has had its telomere repeats removed and each chromosome end trimmed to
    the same fixed length, one record per arm. Producing one is outside this application's
    scope; see the README.
    """
    info = ReferenceInfo()
    records = list(_iter_fasta_lengths(path))
    if not records:
        info.errors.append("no FASTA records found")
        return info

    info.n_records = len(records)
    info.contigs = [name for name, _ in records]

    # --- uniform record length -------------------------------------------------------------
    lengths = {length for _, length in records}
    if len(lengths) == 1:
        info.arm_length_bp = lengths.pop()
    else:
        by_length: dict[int, int] = {}
        for _, length in records:
            by_length[length] = by_length.get(length, 0) + 1
        commonest = max(by_length, key=lambda k: by_length[k])
        odd = [f"{n} ({length:,} bp)" for n, length in records if length != commonest]
        info.errors.append(
            f"records are not all the same length: {len(by_length)} distinct lengths, "
            f"most are {commonest:,} bp. Odd ones: {', '.join(odd[:max_report])}"
            + (f" and {len(odd) - max_report} more" if len(odd) > max_report else "")
            + ". Pipeline A subtracts the arm length from q-arm coordinates, so a ragged "
            "reference corrupts q-arm telomere lengths while leaving p arms correct."
        )

    # --- contig naming ---------------------------------------------------------------------
    bad = [name for name in info.contigs if not CONTIG_RE.match(name)]
    if bad:
        hint = ""
        if any(re.match(r"^chr([1-9])([pq])", n) for n in bad):
            hint = (" Autosomes must be zero-padded (chr1p -> chr01p): the pipeline's own "
                    "padding fix-up is dead code, so unpadded names sort wrongly throughout.")
        info.errors.append(
            f"{len(bad)} of {len(records)} contig names are not of the form "
            f"chr01p / chr22q_MATERNAL: {', '.join(bad[:max_report])}"
            + (f" and {len(bad) - max_report} more" if len(bad) > max_report else "")
            + "." + hint
        )

    # --- arms and haplotypes ---------------------------------------------------------------
    haplotypes = sorted({
        m.group(4) for n in info.contigs if (m := CONTIG_RE.match(n)) and m.group(4)
    })
    info.haplotypes = haplotypes
    if len(haplotypes) > 1:
        info.warnings.append(
            f"reference mixes haplotypes ({', '.join(haplotypes)}). That is supported, but "
            "reads from near-identical chromosome ends may map ambiguously between them."
        )

    arms = {n[:-len(m.group(3))] if (m := CONTIG_RE.match(n)) and m.group(3) else n
            for n in info.contigs if CONTIG_RE.match(n)}
    missing_pairs = sorted(
        a[:-1] for a in arms if a.endswith("p") and a[:-1] + "q" not in arms
    )
    if missing_pairs:
        info.warnings.append(
            f"{len(missing_pairs)} chromosome(s) have a p arm but no q arm: "
            f"{', '.join(missing_pairs[:max_report])}"
        )
    return info


@dataclass
class FastqInfo:
    n_reads: int = 0
    total_bases: int = 0
    gzipped: bool = False
    errors: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def sniff_fastq(path: Path, *, limit: int | None = None) -> FastqInfo:
    """Confirm the file really is FASTQ and count reads.

    `limit` stops after that many reads, for a quick check on a large file.
    """
    info = FastqInfo(gzipped=path.suffix == ".gz")
    try:
        with _open(path) as handle:
            first = handle.readline()
            if not first:
                info.errors.append("file is empty")
                return info
            if not first.startswith("@"):
                info.errors.append(
                    f"does not look like FASTQ: first line starts with {first[:1]!r}, expected '@'"
                )
                return info
            handle.seek(0)
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline()
                plus = handle.readline()
                qual = handle.readline()
                if not qual:
                    info.errors.append(f"truncated record at read {info.n_reads + 1}")
                    break
                if not plus.startswith("+"):
                    info.errors.append(
                        f"malformed record at read {info.n_reads + 1}: third line is not '+'"
                    )
                    break
                if len(seq.strip()) != len(qual.strip()):
                    info.errors.append(
                        f"sequence and quality lengths differ at read {info.n_reads + 1}"
                    )
                    break
                info.n_reads += 1
                info.total_bases += len(seq.strip())
                if limit is not None and info.n_reads >= limit:
                    break
    except (OSError, UnicodeDecodeError, EOFError) as exc:
        info.errors.append(f"could not read: {exc}")
    return info


def guess_barcode(filename: str) -> str | None:
    """The reference corpus names files HG002.1.NB65uq.fastq, so the tag is usually there."""
    match = re.search(r"\b(NB\d{2}uq)\b", filename)
    return match.group(1) if match else None
