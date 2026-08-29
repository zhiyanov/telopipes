# telomers

A small web application for measuring telomere length from Oxford Nanopore reads, for use
within our lab. It wraps two published pipelines from Karimian et al., *Science* 2024:

| Pipeline | Method | Reference genome | Per-chromosome results |
|---|---|---|---|
| `telomere-r` | Mapping-based: `minimap2` → soft clip beyond the subtelomere boundary (R/Bioconductor) | required, pre-cut | yes |
| `telonp` | Sequence-based: per-read boundary detection (GreiderLab TeloBP) | not needed | no |

Each pipeline runs in its own container, orchestrated by Snakemake **inside** that container,
so nothing bioinformatics-related is installed on the host.

## Status

Phases 0–2 are complete: both pipelines run from the command line, reproduce the published
reference results byte-for-byte, and are normalised into a common schema with figures. There
is no web interface yet.

- [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md) — what the host can and cannot do, with evidence
- [`docs/PATCHES.md`](docs/PATCHES.md) — every deviation from the upstream scripts, and how it was verified
- [`docs/PLAN.md`](docs/PLAN.md) — the design and the remaining phases

## Usage

```bash
uv sync
telomers doctor                 # check podman, images, disk
telomers images build           # first build compiles Bioconductor from source; ~5 min
telomers pipelines              # what is available

telomers run telonp \
    --reads sample.fastq.gz --sample-label mysample

telomers run telomere-r \
    --reads sample.fastq.gz --reference cut.MATERNAL.fasta \
    --sample-label mysample --barcode-name NB65uq
```

Registering inputs is optional but hardlinks them, so cataloguing files already on the machine
costs no disk:

```bash
export TELOMERS_IMPORT_DIRS='["/data/nanopore"]'   # registering by path is off until you say where
telomers datasets   add --registered /data/nanopore/sample.fastq
telomers references check  cut.MATERNAL.fasta      # dry run: is this reference usable?
telomers references add --registered cut.MATERNAL.fasta

telomers runs                      # what has been run
telomers results <run-id> --by-arm # canonical summary
telomers reindex                   # rebuild the index by walking the data root
```

Each run gets a self-describing directory under `data/runs/<id>/`:

```
run.json     what produced it: pipeline, image, inputs, parameters
run.log      the container's output
work/        everything the pipeline itself wrote, including its own figures
results/     per_read.csv, summary.json, per_arm.csv and three figures
```

Because every fact also lives on disk, the SQLite index is a cache — `telomers reindex`
rebuilds it. That is why there are no database migrations.

Re-running with the same `--run-id` resumes: Snakemake skips completed steps, so an interrupted
pipeline-A run picks up at the analysis without repeating the alignment.

### A note on the reference genome

Pipeline A needs a *cut* reference: telomere repeats removed, every chromosome end trimmed to
the same length, one record per arm, named `chr01p_MATERNAL`. Producing one is out of scope
here — use TeloBP's `trimGenome.py` — but `telomers references check` will tell you whether
yours is usable and exactly what is wrong if not. This matters more than it looks: a reference
whose records are not all the same length silently corrupts q-arm telomere lengths while
leaving p arms perfect, and unpadded contig names (`chr1p` rather than `chr01p`) do not error
at all, they just sort wrongly throughout.

`scripts/reproduce-reference.sh` re-runs both pipelines on the published HG002 dataset and
diffs every output against `reference/results/`. It is the regression test for any change to
the pipeline scripts.

## Architecture note

Images are built **natively for the host architecture**. On this workstation that is arm64:
the kernel uses 16 KB pages, where x86_64 emulation (FEX, qemu-user) cannot work at all. This
is why pipeline A is based on `rocker/r-ver` rather than the reference's amd64-only
`bioconductor/bioconductor_docker`.

## Building the images

```bash
podman build --platform linux/arm64 -t telomers/pipeline-a:0.1.0 pipelines/telomere-r
podman build --platform linux/arm64 -t telomers/pipeline-b:0.1.0 pipelines/telonp
```

Pipeline A compiles Bioconductor from source on arm64 and takes a while; this is a one-time cost.

## Provenance

Wrapper scripts under `pipelines/*/scripts/` are vendored from the reference implementation and
from the upstream GreiderLab repositories. Every deviation from upstream is recorded in
`docs/PATCHES.md`. TeloBP is installed from GitHub at a pinned commit; see
`pipelines/telonp/requirements.txt`.
