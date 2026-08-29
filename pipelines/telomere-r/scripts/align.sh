#!/usr/bin/env bash
# Method 1, step 2: map ONT reads to the cut (telomere-removed, 500 kb chr-end)
# reference and produce a sorted, indexed BAM. Runs inside the telomere-r image.
#
# Usage: align.sh <reference.fasta> <reads.fastq> <out_bam> [threads]
set -euo pipefail

REF="$1"
FASTQ="$2"
OUT_BAM="$3"
THREADS="${4:-4}"

mkdir -p "$(dirname "$OUT_BAM")"

echo "[align] minimap2 -ax map-ont  (ref=$REF  reads=$FASTQ)"
# -ax map-ont: ONT preset, as specified in the paper's methods.
# --MD keeps mismatch info; -Y uses soft-clipping for supplementary alignments
# so the R softclip-based telomere length logic sees clips consistently.
minimap2 -t "$THREADS" -ax map-ont -Y "$REF" "$FASTQ" \
    | samtools sort -@ "$THREADS" -o "$OUT_BAM" -

samtools index "$OUT_BAM"

echo "[align] done -> $OUT_BAM"
echo "[align] flagstat:"
samtools flagstat "$OUT_BAM"
