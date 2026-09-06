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

## P5 — `run_analysis.R`: named arguments, literal barcodes, stable filenames

Three changes to our own driver, none of which alter results.

**Named arguments** (`optparse`) replace the five positional ones. `readLnMin`, `alnScTh`,
`mapStartTh` and `seqLn` were hardcoded in the call to `alignWriteTableFn()`; they are now
`--read-len-min`, `--aln-sc-th`, `--map-start-th` and `--seq-ln`, which is what lets the
Snakefile and the eventual web form expose them.

**`--barcode-seq` accepts a literal sequence.** Upstream resolved a *name* with `get()`, which
works only for the 18 TeloTags defined in `telobp_functions.R` and otherwise fails with a bare
`Error: object 'NB99uq' not found` from inside R. A name is still accepted via
`--barcode-name`, but an unknown one now produces a message listing the known tags and
pointing at `--barcode-seq`.

**Stable output filenames.** Upstream baked the experiment name, tag name and threshold into
`reads.table_<exp>.prim.Samp.<bc>_alnScTh.20.txt`. The driver renames the produced files to
`reads.table.txt`, `reads.table.all.txt`, `summary.overall.csv` and `summary.by_chr.csv` so
Snakemake can declare them as outputs without globbing. Renaming in the driver rather than
changing `telobp_functions.R` keeps the shared upstream file closer to its original.

**Verified:** reproduction clean, running with `--barcode-seq TTCTCAGTCTTCCTCCAGACAAGG` — so
the literal-sequence path is confirmed to give the same answer as the name lookup.

## P6 — `--emit-all`: write the reads the filters reject

`writeTelTableFn()` computes `indAlnScTh`, `indMapPosTh` and `indPosTel` per read, then writes
out only the rows where all three hold. The rejected reads, and the reason, were discarded —
so read accounting could only be scraped from stdout.

With `--emit-all` the same table is written for every primary read, plus `pass_aln_score`,
`pass_map_start`, `pass_pos_telomere`, `qc_status` and `qc_reason` (the first failing filter,
in the order the pipeline conceptually applies them). On the reference dataset:

```
4732 primary reads (>= 3 kb)
  4362  pass                            <- byte-identical to the tagged table
   347  negative_or_missing_telomere
    21  map_start_too_internal
     2  barcode_below_threshold
```

This is what makes the read funnel a first-class result rather than a log message, and it puts
pipeline A on the same footing as pipeline B, which reports its failures for free via negative
sentinels.

**Verified:** the tagged table is unchanged, and the `pass` rows of the all-reads table match
it exactly (checked by `scripts/reproduce-reference.sh`).

---

## Effect of the memory patches

P2 and P3 together took the R analysis step from **20 s to 7.6 s** on the reference dataset,
and the BAM from **566 MB to 98 MB**, with every output byte-identical.

---

## Reference cutting (new, not a patch to upstream)

`pipelines/telomere-r/scripts/cut_reference.py` derives a cut reference from a whole genome.
Written from scratch rather than calling TeloBP's `trimTeloReferenceGenome()`, which emits
`_C`/`_G` suffixes rather than `p`/`q` and, when an end's telomere length is zero, evaluates
`record[-0-N:-0]` — that is `record[-N:0]` — producing an empty record.

Two decisions worth recording:

**Naming is normalised here, not downstream.** `chr1_MATERNAL` becomes `chr01p_MATERNAL` and
`chr01q_MATERNAL`. This is the only place it can be fixed: `posEndPQchrFn()` zero-pads on one
line and overwrites the result from `rname` on the next, so its own repair never runs, and an
unpadded name does not error — it just sorts wrongly everywhere.

**A candidate telomere tract needs a quality floor.** Repeat matches are walked outward from
the terminus and merged across gaps up to 250 bp, mirroring Telo-seq's
`build_telomere_reference.py`. Without a floor, two chance hexamers near a telomere-less
terminus trim a spurious ~150 bp, which silently shifts that arm's coordinate frame relative to
every other arm. A tract must now be ≥100 bp and ≥60% matched. Five-mers (`CCTAA`, `TTAGG`)
were dropped from the patterns for the same reason: they occur by chance every ~1 kb.

**Validation.** Cutting HG002 v1.2 (maternal) and re-running pipeline A on the same reads
reproduces the v0.7 result almost exactly — median 4,593.5 vs 4,589.5 bp, n 4,358 vs 4,362,
and 36 of 46 arms agree within 50 bp. The two arms that move are both acrocentric: `chr15p`
loses 166 reads and `chr22p` gains 165, i.e. v1.2 separates two near-identical short arms that
v0.7 conflated. `chr14p`, the shortest-telomere result, is unchanged.
