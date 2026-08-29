# Pipeline A (telomere-r): mapping-based telomere length, Karimian et al. Science 2024.
#
# Built NATIVELY for the host architecture. The reference implementation used
# bioconductor/bioconductor_docker, which publishes amd64 only; this host is arm64 with a
# 16 KB page size, where x86_64 emulation (FEX / qemu-user) cannot work. rocker/r-ver is
# multi-arch, so we install the same Bioconductor packages onto it instead.
# See docs/ENVIRONMENT.md.
FROM docker.io/rocker/r-ver:4.4.1

# minimap2 + samtools from the Ubuntu repos, plus python3 for Snakemake.
RUN apt-get update && apt-get install -y --no-install-recommends \
        minimap2 \
        samtools \
        python3 \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

# rocker points R at Posit Package Manager, which serves binaries for amd64 only.
# On arm64 that silently degrades; use CRAN directly so behaviour is the same either way.
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' \
      > /usr/local/lib/R/etc/Rprofile.site

# Bioconductor 3.19 pairs with R 4.4. Biostrings / GenomicAlignments / S4Vectors are the
# three the pipeline actually uses; optparse is for our patched run_analysis.R CLI.
# writexl is deliberately NOT installed: telobp_functions.R loads it but never calls it.
RUN R -q -e "install.packages(c('BiocManager','RColorBrewer','optparse'))" \
 && R -q -e "BiocManager::install(version='3.19', ask=FALSE)" \
 && R -q -e "BiocManager::install(c('Biostrings','GenomicAlignments','S4Vectors'), update=FALSE, ask=FALSE)" \
 && R -q -e "stopifnot(all(c('Biostrings','GenomicAlignments','S4Vectors','RColorBrewer','optparse') %in% rownames(installed.packages())))"

# Snakemake orchestrates this pipeline's rules from inside the container, so the host
# never needs it. Ubuntu 24.04 is PEP 668 managed, hence the venv.
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir snakemake==9.26.1
ENV PATH=/opt/venv/bin:$PATH

COPY scripts/ /opt/telomers/
WORKDIR /work
