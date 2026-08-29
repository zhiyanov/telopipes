# Deviations from the upstream scripts

Every change we make to the vendored GreiderLab code, and why. Each entry states how it was
verified. The rule is that no patch lands until `scripts/reproduce-reference.sh` still
reproduces `reference/results/` byte-for-byte, or until the deviation is deliberate and
explained here.

Upstream sources:

- `pipelines/telomere-r/scripts/telobp_functions.R` — GreiderLab/Telomere-Read-Analysis
  @ `02911abf`, i.e. `human.telomere.pipeline.v0.9.R` with its trailing 25-line exec block
  removed so the file can be `source()`d. The pristine original is kept alongside it at
  `scripts/upstream/human.telomere.pipeline.v0.9.orig.R` for diffing.
- `pipelines/telomere-r/scripts/align.sh`, `run_analysis.R` — written for the reference
  project, not upstream.
- TeloBP — GreiderLab/TeloBP @ `05bd9d8`, installed from git, not vendored. Verified
  byte-identical to the copy that produced `reference/results/` (see `ENVIRONMENT.md`).

---

## P1 — `telobp_functions.R`: drop `library(writexl)`

`writexl` is loaded at the top of the file but `write_xlsx()` is never called anywhere in it;
the per-read table is written with `write.table()`. Removing the load lets us drop the package
from the image.

**Verified:** reproduction clean.

Worth noting the order this had to happen in. The package was dropped from the image first,
while the script was still unpatched, and the run failed at `library(writexl)`. "Loaded but
never called" was accurate, but the `library()` call still needs the package present — so the
image change is only safe together with the source change.

## P2 — `telobp_functions.R`: give `scanBam()` a `ScanBamParam`

Upstream called `scanBam(bamFN, nThreads=)` with no parameters, which materialises every field
of every record into R. `qual` alone is the same size as `seq`, and grep confirms `qual`,
`mrnm`, `mpos` and `isize` appear **zero** times in the file — only `qname`, `flag`, `rname`,
`strand`, `pos`, `qwidth`, `mapq`, `cigar` and `seq` are used. `seq` **is** required:
`align2tagPQchrFn()` scores the TeloTag barcode against the read ends.

**Verified:** reproduction clean.

## P3 — `align.sh`: drop secondary and supplementary records before sorting

`minimap2 -Y` forces soft clipping on supplementary alignments, so each carries a full copy of
`SEQ`. On the reference dataset that is most of the file:

```
33372 total records   5961 primary   17147 secondary   10264 supplementary
```

82% of records were secondary or supplementary, and every one of them was being written to
disk and then loaded into R. Filtering with `samtools view -u -F 0x900` before the sort takes
the BAM from **566 MB to 98 MB (5.8x)** and removes the same volume from `scanBam`'s working
set, compounding with P2.

Nothing downstream consumes them: `processBamFn()` indexes `mapStats$primary` throughout, and
`mapStats$secondary` / `$suppplemental` feed only unused locals. Unmapped reads are **not**
dropped — flag `0x4` is not part of `0x900` — so the read-length histogram is unaffected.
Full flag accounting is preserved by teeing the unfiltered stream through `samtools flagstat`
into `flagstat.txt`.

**Verified:** reproduction clean. `scripts/reproduce-reference.sh` compares
`samtools view -F 0x900` on both sides, which is the like-for-like comparison now that the
BAMs differ by construction.

One caveat, recorded because it is invisible: `reads$chr` is a factor whose levels come from
the rnames actually observed, and the by-chromosome plots iterate those levels. A contig that
received *only* secondary alignments would lose an empty axis category. With a 46-contig
reference and 5,961 reads this cannot realistically arise, and it did not here.

## P4 — `run_telonp.py`: stream and chunk instead of pandas + pandarallel

The reference driver loaded the entire FASTQ into a DataFrame and dispatched with
`pandarallel.parallel_apply`, which pickles a slice per worker **through `/dev/shm`**. Two
problems:

1. Peak memory scaled as *(total sequence bytes x workers)*.
2. podman caps `/dev/shm` at 64 MB, so the reference driver dies with
   `OSError: [Errno 28] No space left on device` in a container — on a host with 83 GB free —
   unless `--shm-size` is raised.

v2 streams reads and scores them in bounded chunks with `ProcessPoolExecutor`, writing CSV rows
as it goes. Peak memory is now *O(chunk_reads x workers)* and no shared memory is used, so no
`--shm-size` is needed. pandas is gone from the compute path entirely, so its version no longer
affects anything (the reference environment had resolved pandas 3.0.5 for a package developed
against 2.2).

Also in v2: `.fastq.gz` is accepted, since `SeqIO.parse(path, "fastq")` does not decompress;
and the `sys.path` surgery is removed. That existed because the vendored TeloBP checkout sat
next to the driver and shadowed the installed package as a namespace package — impossible in
the image, where TeloBP is in site-packages and the driver is in `/opt/telomers`. A failed
import now reports itself instead of silently half-working.

**Verified:** all three outputs byte-identical to `reference/results/telonp/`.
