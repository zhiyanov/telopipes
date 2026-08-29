# Pipeline B (TeloNP): reference-free per-read telomere boundary detection.
# python:3.12-slim is multi-arch, so this builds natively on arm64.
FROM docker.io/library/python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY scripts/ /opt/telomers/
ENV MPLBACKEND=Agg
WORKDIR /work
