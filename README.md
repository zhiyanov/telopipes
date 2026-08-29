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

Early. Phase 0 (environment proof) is complete — see [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md).
The design is in [`docs/PLAN.md`](docs/PLAN.md).

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
