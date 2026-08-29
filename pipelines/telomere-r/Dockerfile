# Pipeline A (telomere-r): mapping-based telomere length, Karimian et al. Science 2024.
#
# Built NATIVELY for the host architecture. The reference implementation used
# bioconductor/bioconductor_docker, which publishes amd64 only; this host is arm64 with a
# 16 KB page size, where x86_64 emulation (FEX / qemu-user) cannot work. rocker/r-ver is
# multi-arch, so we install the same Bioconductor packages onto it instead.
# See docs/ENVIRONMENT.md.
FROM docker.io/rocker/r-ver:4.4.1

# minimap2 + samtools from the Ubuntu repos, plus python3 for Snakemake.
#
# The -dev headers matter here: there are no upstream Bioconductor binaries for arm64, so
# every package compiles from source. Without them the build fails deep in the dependency
# tree with misleading errors -- zlib/bzip2/lzma are needed by XVector and Rhtslib,
# curl+ssl by httr -> UCSC.utils -> GenomeInfoDb, and pcre2/deflate to relink rocker's
# own littler when BiocManager refreshes it.
RUN apt-get update && apt-get install -y --no-install-recommends \
        minimap2 \
        samtools \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libpcre2-dev \
        libdeflate-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# rocker points R at Posit Package Manager, which serves binaries for amd64 only.
# On arm64 that silently degrades; use CRAN directly so behaviour is the same either way.
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' \
      >> /usr/local/lib/R/etc/Rprofile.site

# Bioconductor 3.19 pairs with R 4.4. Biostrings / GenomicAlignments / S4Vectors are the
# three the pipeline actually uses; optparse is for our patched run_analysis.R CLI.
RUN R -q -e "install.packages(c('BiocManager','RColorBrewer','optparse'))" \
 && R -q -e "BiocManager::install(version='3.19', ask=FALSE, update=FALSE)" \
 && R -q -e "BiocManager::install(c('Biostrings','GenomicAlignments','S4Vectors'), update=FALSE, ask=FALSE)" \
 && R -q -e "stopifnot(all(c('Biostrings','GenomicAlignments','S4Vectors','RColorBrewer','optparse') %in% rownames(installed.packages())))"

# Bioconductor 3.19 split pairwiseAlignment() out of Biostrings into its own package. The
# reference ran on Bioc 3.18, where it still lived in Biostrings, so nothing in the upstream
# script asks for pwalign -- but alignWriteTableFn() dies without it when scoring the TeloTag
# barcode. Same code, new home; whether it also behaves identically is what the reproduction
# in scripts/reproduce-reference.sh checks.
RUN R -q -e "BiocManager::install('pwalign', update=FALSE, ask=FALSE); stopifnot('pwalign' %in% rownames(installed.packages()))"

# Snakemake orchestrates this pipeline's rules from inside the container, so the host never
# needs it. rocker ships no Python, and the base's Ubuntu release would dictate which version
# apt gives (r-ver:4.4.x is 22.04 -> Python 3.10, too old for Snakemake 9). Install a
# standalone interpreter with uv instead, so the orchestrator does not depend on whichever
# rocker base we happen to sit on -- and stays on the same Snakemake as pipeline B.
COPY --from=ghcr.io/astral-sh/uv:0.9.9 /uv /usr/local/bin/uv
ENV UV_PYTHON_INSTALL_DIR=/opt/uv-python
RUN uv venv --python 3.12 /opt/venv \
 && VIRTUAL_ENV=/opt/venv uv pip install --no-cache snakemake==9.26.1
ENV PATH=/opt/venv/bin:$PATH

COPY scripts/ /opt/telomers/
# Snakemake runs inside this image, so its workflow ships in it too.
COPY Snakefile /opt/telomers/Snakefile
WORKDIR /work
