# Telo-seq analysis with chromosome assignment.
#
# FROM the reference-free image, so NCRF, its converted helpers and the length caller are
# shared rather than duplicated -- the same relationship the upstream project uses
# (Dockerfile.teloseq is FROM telomere-r). Build pipeline-c first; `telomers images build`
# does so because "teloseq" sorts before "teloseq-mapped".
FROM localhost/telomers/pipeline-c:0.1.0

RUN apt-get update && apt-get install -y --no-install-recommends \
        minimap2 \
        samtools \
    && rm -rf /var/lib/apt/lists/*

COPY scripts/ /opt/telomers/
COPY Snakefile /opt/telomers/Snakefile
WORKDIR /work
