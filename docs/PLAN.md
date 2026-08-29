# Telomere Analysis Web App — Design Plan (simplified)

## Context

`/home/antzhi/bio/reference` reproduces Karimian et al. (*Science* 2024) telomere-length
measurement on HG002 ONT reads, but as a flat pile of scripts driven by one `run_all.sh` with
every path, sample, barcode and parameter hardcoded. Analysing a new sample means editing the
shell script.

This repo turns the two pipelines that matter into a small web app for **the wet-lab members of
our own group**: register reads and a reference, pick a pipeline, run it, browse results.

- **Pipeline A — `telomere-r`** (mapping-based): `minimap2 -ax map-ont -Y` against a *cut*
  reference; telomere = soft-clipped bases beyond the subtelomere boundary. R + Bioconductor.
  Two steps. Gives per-chromosome-arm results.
- **Pipeline B — `TeloNP`** (sequence-based): reference-free per-read boundary detection via
  GreiderLab TeloBP. Pure Python, one step. Bulk distribution only.

Deferred: Telogator2, Telo-seq, cross-pipeline comparison.

## What was cut from the first draft, and why

A handful of people on a LAN, one run at a time, on one machine. That removes most of the
infrastructure the first draft carried:

| Dropped | Replaced by | Why it was safe |
|---|---|---|
| Auto-cutting a full genome into arm form | Upload an already-cut reference | Your call — removes the trickiest, most bug-prone component (see the `trimTeloReferenceGenome` empty-record bug) and a whole prep-job code path |
| Custom `Runner` protocol + 3 implementations + `--detach` + crash recovery + `reconcile`/`adopt` | Snakemake, running inside the pipeline container | Your call. The app makes one `podman run` per run. If it dies mid-run the run is marked interrupted and **Resume** re-issues the same command; Snakemake skips completed outputs |
| Bespoke step cache (`stepcache/`, cache keys, hardlinking) | Snakemake's incremental execution | This is exactly what Snakemake already does — it was the strongest argument for using it |
| Content-addressed blob store, dedup, refcounting | One directory per dataset/reference | Nobody is uploading the same 250 MB FASTQ twice on purpose; sha256 is kept for integrity display only |
| tus-lite resumable uploads + ~60 lines of JS | Plain streaming upload + register-by-path | LAN upload of a few hundred MB; register-by-path is the path that will actually get used |
| Alembic migrations, `artifact` table | `create_all`, artifacts discovered by scanning the run dir | The DB becomes a **rebuildable index** — every fact lives in `run.json`/`meta.json` on disk, so `telomers reindex` recreates it. No migration story needed |
| Spec-driven Pydantic model generation, `render.py`, declarative `tabular` adapter | Two hand-written ~30-line params models in a registry dict | With two pipelines, generating models from a spec is more code than writing them |
| SSE log streaming, `--log-driver k8s-file` | HTMX polling a log file by byte offset | One update mechanism instead of two |
| Plot cache, filter UI, SVG pipeline | Three PNGs rendered once after the run | Filters can come later if anyone asks |
| Second "native arm64" image variant per pipeline | One image per pipeline | Moot after phase 0 — everything is native arm64 now; there is no amd64 variant to fall back from |
| Retention policies, disk widgets, pre-flight refusal, quadlet, auth hooks, `--network=none`/`--pids-limit`/tmpfs tuning | A "delete run" button and a `--memory` cap | Two users, 85 GB free |

Kept because they protect against silent wrong numbers, not because they're infrastructure:
the **reference validator**, the **golden-file tests against the reference outputs**, and the
**R memory fixes**.

## Decisions

| | |
|---|---|
| Job execution | **Snakemake**, one workflow per pipeline, running **inside** that pipeline's container |
| Containers | podman, **native arm64 images** — x86 emulation is impossible here (see `docs/ENVIRONMENT.md`) |
| Stack | FastAPI + Jinja2 + HTMX, server-rendered, no node |
| Reference genome | **User supplies an already-cut reference**; we validate it, never modify it |
| Barcode | Dropdown of the 18 known TeloTags plus a custom sequence |
| Storage | Plain directories + SQLite index that can be rebuilt from disk |
| Host footprint | 38-package venv in the project dir; nothing system-wide |
| Vendoring | Vendor our wrapper scripts; TeloBP from GitHub at a pinned commit |

### Environment (probed)

podman 5.8.4 rootless, overlay driver; no docker/compose/node; uv 0.12.7; 8 cores / 15 GB /
~85 GB free; btrfs; **SELinux Enforcing**; `podman.socket` present but disabled (we shell out
to the CLI, so it stays that way).

**Phase 0 killed the amd64-under-emulation plan.** The kernel is
`7.1.6-400.asahi.fc44.aarch64+16k` with a **16 KB page size**. x86_64 guarantees 4 KB pages, so
FEX crashes with `<jemalloc>: Unsupported system page size`, and qemu-user hits the same wall —
it is a kernel property, not a misconfiguration. Native arm64 was verified working and every
dependency is available for it (`minimap2`, `samtools` in the distro repos; `rocker/r-ver` and
`python:3.12-slim` multi-arch). Only `bioconductor/bioconductor_docker` is amd64-only, so
pipeline A changes base image. Full evidence in `docs/ENVIRONMENT.md`.

---

## 1. Layout

```
telomers/
├── Makefile  pyproject.toml  uv.lock  .env.example
├── docs/  ENVIRONMENT.md   PATCHES.md   (deviations from upstream, with rationale)
├── pipelines/
│   ├── telomere-r/
│   │   ├── Snakefile  Dockerfile  meta.toml
│   │   └── scripts/  align.sh  run_analysis.R  telobp_functions.R
│   │                 upstream/human.telomere.pipeline.v0.9.orig.R   (pristine, for diffing)
│   └── telonp/
│       ├── Snakefile  Dockerfile  requirements.txt
│       └── scripts/  run_telonp.py
├── src/telomers/
│   ├── cli.py  config.py  storage.py  db.py  doctor.py
│   ├── pipelines.py        registry: id → {meta, snakefile, params model, normalizer}
│   ├── params.py           two Pydantic models (TelomereRParams, TeloNPParams)
│   ├── execute.py          write config.yaml, one `podman run`, capture log, cancel
│   ├── validate.py         FASTQ sniffing + cut-reference validator
│   ├── results.py          normalize → per_read.csv / summary.json / per_arm.csv
│   ├── plots.py            three matplotlib PNGs
│   └── web/  app.py  routes.py  templates/  static/vendor/htmx.min.js
└── tests/  unit/  data/{golden/, refs/, tiny.fastq}
```

`src/telomers/` is a src-layout package; console script `telomers`.

`telobp_functions.R` is **vendored**, not fetched at build time: it is upstream
`human.telomere.pipeline.v0.9.R` with its trailing exec block stripped, and `source()`ing the
pristine file would execute that block. Provenance: `GreiderLab/Telomere-Read-Analysis` @
`02911abf`.

---

## 2. Pipelines as Snakemake workflows

One `Snakefile` per pipeline, **baked into that pipeline's image and executed inside it**.

Every method in and out of scope happens to be a *single-image* pipeline: `telomere-r` carries
minimap2, samtools and R together, so both of its steps run in one container; TeloNP is one step;
and the deferred Telogator2 and Telo-seq images are self-contained too (the latter is literally
`FROM telomere-r`). Nothing needs an orchestrator that spans images.

That means Snakemake never touches the host. The app makes one plain `podman run` per run and
Snakemake handles the DAG, incremental re-execution and logging from inside. No podman socket, no
nested podman, no host/container path mismatch — the design has no component reaching across the
container boundary.

```python
# pipelines/telomere-r/Snakefile   → COPYed to /opt/telomers/Snakefile
configfile: "/work/config.yaml"

rule all:
    input: "r_analysis/reads.table.txt", "r_analysis/summary.by_chr.csv"

rule align:
    input:  reads="/inputs/reads.fastq", ref="/refs/reference.fasta"
    output: bam="aligned.sort.bam", bai="aligned.sort.bam.bai", flagstat="flagstat.txt"
    threads: config["threads"]
    shell:  "bash /opt/telomers/align.sh {input.reads} {input.ref} {output.bam} {threads}"

rule analyze:
    input:  bam="aligned.sort.bam"
    output: "r_analysis/reads.table.txt", "r_analysis/reads.table.all.txt",
            "r_analysis/summary.overall.csv", "r_analysis/summary.by_chr.csv"
    params: cfg=config
    shell:
        "Rscript /opt/telomers/run_analysis.R --bam {input.bam} --outdir r_analysis"
        " --exp-name '{params.cfg[sample_label]}' --chr-arm-ln {params.cfg[chr_arm_ln]}"
        " --barcode-seq {params.cfg[barcode_seq]} --read-len-min {params.cfg[read_len_min]}"
        " --aln-sc-th {params.cfg[aln_score_threshold]}"
        " --map-start-th {params.cfg[map_start_threshold]} --emit-all"
```

Plain shell commands — no `podman run` string-building, and no fighting Snakemake's `{}`
formatting. `pipelines/telonp/Snakefile` is one rule with no reference input.

Output filenames are **fixed**, not the reference's
`reads.table_<exp>.prim.Samp.<bc>_alnScTh.20.txt`, so rules and the normalizer never glob. That
is part of the patch set below.

`meta.toml` holds display metadata only — title, citation,
`provides = ["per_read","bulk","per_arm"]`, `requires_reference = true`, image tag. Parameters are
a hand-written Pydantic model in `params.py`; the run form is rendered from its JSON schema.

### Patches to the upstream scripts (all recorded in `docs/PATCHES.md`)

1. **`run_analysis.R`: `optparse` named args**, and take `--barcode-seq` (a literal sequence)
   instead of resolving a name with `get(barcodeName)`. This is what makes the custom-barcode
   option possible — today an unknown name gives a bare `Error: object 'NB99uq' not found`.
2. **`run_analysis.R`: stable output filenames**, so Snakemake can declare them.
3. **`run_analysis.R` / `writeTelTableFn`: `--emit-all`**, writing `reads.table.all.txt`
   including failed reads with the per-read filter flags the code already computes
   (`indAlnScTh`, `indMapPosTh`, `indPosTel`). ~6 lines, and it turns the read funnel
   (5,961 → 4,929 → 4,362) into real data instead of a stdout scrape.
4. **`align.sh`: `samtools view -b -F 0x900`** before sorting, with flagstat preserved via a
   `tee`. `minimap2 -Y` forces soft-clipping on supplementary alignments, so they carry full
   SEQ — which is why 5,961 reads produce a 593 MB BAM. Prefiltering shrinks it ~10–20×.
5. **`telobp_functions.R` / `processBamFn`: `ScanBamParam(what=...)`.** `scanBam` is currently
   called with no params, loading every field of every record. Verified by grep: `qual`,
   `mrnm`, `mpos`, `isize` appear **zero times** in the file, and `qual` alone is the same byte
   count as `seq`. (`seq` **is** required — the barcode alignment uses it.)
6. **`run_telonp.py` v2:** `gzip.open` for `.fastq.gz`; stream reads into a
   `ProcessPoolExecutor` in chunks and write CSV incrementally, dropping pandas and pandarallel
   from the compute path. Peak RSS goes from `sequence_bytes × (1 + workers)` — the reference's
   fork-and-pickle model, and pipeline B's dominant limit — to `O(chunk × workers)`. Keep
   `matplotlib.use("Agg")`; drop the `sys.path` surgery (in the image TeloBP is in
   site-packages and the driver is in `/opt/telomers/`, so the shadowing it guarded against
   can't occur).

Patches 4 and 5 are output-neutral and get verified **once** in phase 0 by diffing
`reads.table` against the reference's committed output — not by a permanent test matrix or a
runtime toggle.

---

## 3. Job execution

The whole layer is one module, `src/telomers/execute.py`, ~120 lines. Per run it writes
`config.yaml` and makes a single call:

```
podman run --rm --platform linux/arm64 --pull=never
  --name telomers-<run_id> --memory <cap> --cpus <n>
  -v <run>/work:/work:z
  -v <dataset>/reads.fastq:/inputs/reads.fastq:ro,z
  -v <reference>/reference.fasta:/refs/reference.fasta:ro,z
  -w /work  telomers/pipeline-a
  snakemake --snakefile /opt/telomers/Snakefile --cores <n> --rerun-incomplete
```

- Launched with `asyncio.create_subprocess_exec`, stdout+stderr teed into `<run>/run.log`, which
  the web UI tails.
- **Concurrency: one run at a time** (`asyncio.Semaphore(1)`). Pipeline A's analyze step budgets
  ~10 of 15 GB and nobody in the lab needs two at once. Queued runs sit in the DB.
- **Cancel:** `podman stop --time 10 telomers-<run_id>`, then `podman kill`. One container, one
  name — no process trees or label sweeps.
- **Resume:** re-issue the identical command. Snakemake's `.snakemake/` state directory lives in
  `/work`, which is the bind-mounted run directory, so it survives the container exiting; an
  interrupted pipeline-A run resumes at `analyze` without redoing the alignment. On startup any
  run left `running` is marked `interrupted` and offered a **Resume** button.
- **Failure:** non-zero exit → run `failed`, `error_summary` = last 40 log lines. Partial outputs
  stay in `work/` and remain browsable. An OOM kill is reported specifically ("the R analysis step
  exceeded its memory cap") because it is the expected failure mode for pipeline A on a large
  sample.
- Every bind mount carries **`,z`** — SELinux is Enforcing here, and a missing label shows up as
  `Permission denied` *inside* the container, which reads like a pipeline bug.
- We mount only the input files read-only and the run's own `work/` read-write. Never the repo,
  never the data root — the reference's `-v $ROOT:/work` hands the container write access to every
  input and every prior result.

**Host footprint:** the app venv resolves to **38 packages** (FastAPI, uvicorn, Jinja2, SQLModel,
pydantic-settings, matplotlib, pandas and their transitives). Snakemake's 40-package tree lives in
the images, along with everything bioinformatics — R, Bioconductor, minimap2, samtools, TeloBP.
Nothing is installed system-wide; `dnf list installed` never changes.
---

## 4. Storage and index

```
$TELOMERS_DATA_ROOT/
├── telomers.sqlite3                        rebuildable index
├── datasets/<id>/    reads.fastq[.gz]  meta.json
├── references/<id>/  reference.fasta   meta.json      (arm_length_bp, contigs, n_records)
└── runs/<id>/
    ├── run.json          pipeline, params, input ids, status, timings
    ├── config.yaml       what Snakemake was given
    ├── snakemake.log
    ├── work/             everything the pipeline wrote (BAM, R tables, its own PNGs)
    └── results/          per_read.csv, summary.json, per_arm.csv, three PNGs
```

**SQLite via SQLModel, `create_all`, no migrations** — because every fact also lives on disk in
`run.json`/`meta.json`, the database is a cache. `telomers reindex` rebuilds it by walking the
data root. That property is what lets us skip Alembic without risking the run history.

Three tables: `dataset`, `reference`, `run`. Artifacts are **not** a table — the results page
lists `work/` and `results/` directly. `run.inputs` is a small JSON column
(`{"reads": "<id>", "reference": "<id>"}`), so a future pipeline with different inputs needs no
schema change.

**Ingest, two ways:**
1. **Register by path** — the one that will actually be used. Point at a file under a
   `TELOMERS_IMPORT_DIRS` allowlist (e.g. `/home/antzhi/bio/reference/data`), hardlink it in when
   on the same filesystem. Zero copy for the existing 3.1 GB corpus. The allowlist matters: this
   is a local-file-read primitive.
2. **Browser upload** — stream `request.stream()` straight to disk with an incremental sha256.
   Never buffer a 500 MB body, and don't use `UploadFile` (it double-writes through a spool file).

`.fastq.gz` accepted; minimap2 reads it natively and `run_telonp.py` v2 handles it.

### Reference validation (kept — this is the main footgun guard)

The user supplies an already-cut reference. We validate and **never modify** it:

- **Every record exactly the same length**, and `arm_length_bp` is **derived from the file**, never
  entered by the user. q-arm coordinates are computed as `pos - chrArmLn - 2`, so a mismatch
  corrupts **q-arm lengths only** while p arms stay perfect — a half-broken failure that looks
  like biology.
- **Contig names match** `^chr(0[1-9]|1[0-9]|2[0-2]|X|Y)[pq](_(MATERNAL|PATERNAL))?$`. Names are
  load-bearing and the R's own repair is dead code: `posEndPQchrFn:207` zero-pads `chr1p`→`chr01p`,
  then `:208` immediately re-derives `reads$chr` from `reads$rname` and discards it. Unpadded
  names don't error — `grepl("chr([0-9]*|X|Y)p")` still matches — you just silently get a
  wrongly-sorted arm label set.
- Rejections show a per-record diff of what's wrong, so a biologist can fix the FASTA headers.

Because `arm_length_bp` comes from the file, `chr_arm_ln` is displayed read-only on the run form
rather than being a parameter anyone can get wrong.

---

## 5. Results

Per pipeline, a ~40-line normalizer writes into `runs/<id>/results/`.

**`per_read.csv`** — CSV, not parquet; 5–6 k rows doesn't justify pyarrow.

`read_id, telomere_bp, read_length_bp, chrom_arm, haplotype, strand, mapq, qc_status, qc_reason`

Pipeline B leaves `chrom_arm`/`haplotype`/`strand`/`mapq` empty. Failed reads are included with
`qc_status='fail'` — B gives them free via its negative sentinels
(`-1 init`, `-10 fused_read`, `-20 strand_type`, `-1000 seq_not_found`), A via the `--emit-all`
patch.

**No adapter correction, in either pipeline.** Both already measure the same quantity — bases
from the read terminus to the boundary, **including** the ligated TeloTag. Pipeline A's column is
literally `Telomere+Adap`; pipeline B's `getTeloNPBoundary` returns `boundaryPoint = x * windowStep`
from a window loop that starts at index 0 (`seq[i:i+teloWindow]`) or at the last base
(`seq[len(seq)-i-teloWindow : len(seq)-i]`) — verified at `TeloBP.py:86-90`, `:173` — so the tag
is counted, not detected. Subtracting 24 bp from one and not the other would *introduce* a bias.

**`summary.json`**: read funnel, n/mean/median/sd/min/max, `qc_breakdown`, resolved params.
**`per_arm.csv`** (pipeline A only): `chrom_arm, haplotype, n, mean, median, sd`, recomputed by us
from `per_read.csv` rather than parsed from the R's `summary.by_chr.csv` — the R collapses
haplotypes, and recomputing means a mapq filter can be applied later without re-running anything.

`plots.py` renders three PNGs once: bulk histogram with median line; per-arm boxplot ordered
`chr01p`…`chrXq` (pipeline A); a QC funnel bar. The pipeline's own PNGs (~23 from R, 1 from
TeloNP) are kept and listed on the results page as "figures produced by the pipeline" — they're
the citable ones, directly comparable to the paper. Note none is guaranteed: every `png(...)`
block in the R sits inside `try({})`, so the page reports "N figures produced" rather than
assuming 23.

---

## 6. Web

```
GET  /                          runs list + "new run"
     /datasets  /references     list, detail, add (upload or register-by-path)
     /runs/new                  pipeline picker → form from the Pydantic JSON schema
POST /runs                      validate → create → enqueue
GET  /runs/{id}                 status + log tail, HTMX poll every 2 s
     /runs/{id}/results         summary table, three PNGs, pipeline figure gallery
     /runs/{id}/files/{path}    download anything under work/ or results/
POST /runs/{id}/cancel  /resume  /delete
GET  /healthz                   podman + FEX probe
```

**One update mechanism: HTMX polling.** `hx-trigger="every 2s"` on the status fragment; the
server stops emitting the polling attribute once the run is terminal, so it self-terminates with
no client logic. The log tail is the same fragment reading the last N KB of `snakemake.log`.

Tabs are driven by `meta.provides` — pipeline B's per-arm tab is **absent, not empty**, with a
one-line explanation. An empty per-arm chart reads as "this sample has no telomeres."

---

## 7. Images

Each image carries its pipeline's scripts, its `Snakefile`, **and Snakemake itself**, so the
image is the complete, self-describing unit of execution.

**`telomers/pipeline-a`** — `FROM docker.io/rocker/r-ver:4.4.1` (multi-arch; the reference's
`bioconductor/bioconductor_docker` is amd64-only and unusable here). apt `minimap2 samtools
python3 python3-venv`, then `BiocManager::install(version='3.19')` and the same three packages
`Biostrings`, `GenomicAlignments`, `S4Vectors`, plus CRAN `RColorBrewer` and `optparse`.
**Drop `writexl`** — the R loads it and never calls it. Snakemake goes in a venv at `/opt/venv`
(Ubuntu 24.04 is PEP 668 managed). `COPY scripts/ /opt/telomers/`.

Two deviations from the reference worth tracking: **R 4.3 → 4.4 and Bioconductor 3.18 → 3.19**,
and rocker's default Posit Package Manager repo is overridden to CRAN because PPM serves
binaries for amd64 only — on arm64 it would silently fall back to source anyway, so the
Dockerfile makes that explicit. Both are why the smoke test against the published numbers
matters.

**`telomers/pipeline-b`** — `FROM docker.io/library/python:3.12-slim`, pinned
`requirements.txt` (TeloBP deps plus Snakemake), TeloBP from
`git+https://github.com/GreiderLab/TeloBP@05bd9d8bb27b9d52e69e9c1c4b6607744365ea9f` (verified
reachable). Pin `pandas`/`numpy`/`biopython` ourselves: TeloBP's `setup.py` has its constraints
commented out, which is exactly how the reference ended up on pandas 3.0.5 for a package
developed against 2.2. `biopython==1.87` matches what produced `reference/results/`; note
1.83 (TeloBP's era) publishes no arm64 wheel and the slim base ships no compiler.

Scripts and Snakefiles go in by `COPY`, never bind-mounted, so the image identifies what ran.
`make images` → `telomers images build`, a thin wrapper over `podman build --platform linux/amd64`.
Build once: the first build pulls ~2–3 GB of base layers under emulation.

---

## 8. Config

`pydantic-settings`, prefix `TELOMERS_`, `.env`:
`DATA_ROOT`, `DATABASE_URL`, `IMPORT_DIRS` (allowlist), `SNAKEMAKE_CORES`, `STEP_MEMORY`,
`IMAGE_TAG`, `LOG_LEVEL`. That's the whole surface.

---

## 9. Testing

- **Golden-file tests, no container needed.** `reference/results/r_analysis/reads.table_*.txt`
  (417 KB) and `reference/results/telonp/*.telonp.all.csv` (~250 KB) get committed under
  `tests/data/golden/`. The normalizers run against them and assert the real numbers in
  milliseconds on every commit: **n 4362, median 4589.5, mean 4827.2** for A; **n 5049, median
  4378.0, pass_rate 0.8470, err_init 674, err_strandType 228, err_fusedRead 10** for B.
- **Unit:** cut-reference validator against `tests/data/refs/{ok,unpadded,ragged}.fasta`; params
  models (unknown barcode, out-of-range threads); FASTQ sniffing; `config.yaml` generation.
- **A tiny fixture** (~40 reads, mini 50 kb reference) for a manual end-to-end check. Cut the
  mini reference **p arms from the start of the record, q arms from the end** — because
  `Telomere = SoftClipq − |posEnd| + 1` with `posEnd` derived by subtracting `chrArmLn`, cutting
  from the correct end makes tiny-reference lengths identical to the full-reference run.
- **One manual smoke** on `HG002.1.NB65uq`, reproducing the medians above. Not automated in CI.

---

## 10. Remaining risks

1. **Toolchain drift from the reference.** Going native moved pipeline A to R 4.4 /
   Bioconductor 3.19, and distro `minimap2`/`samtools` versions differ from those that produced
   `reference/results/`. Separately, the TeloTag step scores with `Biostrings::pairwiseAlignment`
   against a fractional threshold of 20, so a read sitting exactly on the boundary could classify
   differently. Both are detected by reproducing n=4362 / median=4589.5 — which is precisely what
   the smoke test is for. This replaced the old FEX-performance risk: emulation is off the table,
   and native execution is strictly faster.
2. **R memory.** Patches 4 + 5 above should take the 593 MB BAM down by an order of magnitude.
   The `--memory` cap means a bigger sample fails cleanly instead of OOM-killing the workstation.
3. **The reference must be pre-cut correctly.** Now entirely on the validator, since we no longer
   cut anything. Worth a short "how to prepare a reference" section in the README pointing at
   TeloBP's `trimGenome.py` and the naming rules.
4. **Snakemake-in-the-image assumes single-image pipelines.** True for all four methods here
   (`telomere-r` bundles minimap2, samtools and R; the deferred Telo-seq image is `FROM
   telomere-r`), but a future pipeline needing two *different* images in one workflow couldn't use
   this shape and would need a host-side launcher. Worth re-checking before adding one.
5. **First-run cold start.** Pipeline A's image compiles Bioconductor from source on arm64
   (no upstream binaries for this architecture), so the first build is long. One-time cost;
   `telomers images build` exists so it happens once, deliberately.

---

## 11. Phases

**Phase 0 — Environment proof. COMPLETE — see `docs/ENVIRONMENT.md`.** Established that amd64
emulation cannot work on a 16 KB-page kernel and that native arm64 does; built both images
(5 min / 30 s); confirmed the vendored TeloBP is byte-identical to the pinned upstream commit;
and **reproduced both pipelines byte-for-byte** against `reference/results/` using the unpatched
scripts — pipeline A `n=4362, median=4589.5` with identical alignment records and per-read
table, pipeline B `n=5049, median=4378.0` with all three outputs identical. Two things it
turned up that the design must carry: Bioconductor 3.19 needs `pwalign` for
`pairwiseAlignment()`, and pandarallel needs `--shm-size` because podman caps `/dev/shm` at
64 MB. Superseded phase-0 text: build both
images and time it. Check the Bioconductor image's default USER (rootless write access to `/work`)
and whether it ships `python3`/`pip` for Snakemake.
Diff vendored TeloBP against `05bd9d8`. Re-run the reference's two commands under podman with the
new mount layout and confirm median 4,589.5 reproduces. Apply patches 4 + 5 and diff `reads.table`
against the committed reference output. Measure pipeline B under FEX vs native.
*Deliverable: `docs/ENVIRONMENT.md` with real numbers and a go/no-go.*

**Phase 1 — Pipelines runnable from the CLI.** Patched scripts, both Snakefiles, both images,
`config.yaml` generation, `telomers run <pipeline> --reads … --reference …`, `telomers doctor`.
*This already replaces `run_all.sh`.*

**Phase 2 — Storage, validation, normalization.** Data root layout, SQLite index + `reindex`,
register-by-path, the reference validator, both normalizers, `plots.py`, golden tests.
*`telomers results show <run>` prints the canonical summary.*

**Phase 3 — Web, read-only.** Runs/datasets/references lists and details, results page, figure
gallery, downloads. Browses runs created by the CLI.

**Phase 4 — Web, write path.** Upload, run form from the params schema, submit, the one-at-a-time
worker, status/log polling, cancel, resume, delete. *Full story works in a browser.*

**Phase 5 — Polish.** README with the reference-preparation guide, `/healthz`, tidy errors.

---

## 12. Verification

- `telomers doctor` → podman present, `uname -m` = `x86_64` inside an amd64 container, SELinux
  mode, free disk.
- `make images`; `podman images` lists both.
- `telomers run telonp --reads tests/data/tiny.fastq` finishes in seconds.
- `pytest` green, including the golden-file assertions of the reference medians.
- Manual smoke on `HG002.1.NB65uq` reproduces n=4362/median=4589.5 and n=5049/median=4378.0.
- `make dev` → `localhost:8000`: register a FASTQ and a cut reference by path, run both pipelines,
  watch the log, view the per-arm boxplot, download `per_read.csv`. Kill the server mid-run,
  restart, hit **Resume**, and confirm Snakemake skips the alignment.
