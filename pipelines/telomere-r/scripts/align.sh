#!/usr/bin/env bash
# Method 1, step 1: map ONT reads to the cut (telomere-removed, 500 kb chr-end) reference and
# produce a sorted, indexed BAM.
#
# Usage: align.sh <reference.fasta> <reads.fastq> <out_bam> [threads]
set -euo pipefail

REF="$1"
FASTQ="$2"
OUT_BAM="$3"
THREADS="${4:-4}"

OUT_DIR="$(dirname "$OUT_BAM")"
FLAGSTAT="$OUT_DIR/flagstat.txt"
mkdir -p "$OUT_DIR"

echo "[align] minimap2 -ax map-ont  (ref=$REF  reads=$FASTQ)"
# -ax map-ont: ONT preset, as specified in the paper's methods.
# -Y uses soft-clipping for supplementary alignments so the R softclip-based telomere length
# logic sees clips consistently.
#
# PATCH: drop secondary and supplementary records (-F 0x900) before sorting. Because -Y forces
# soft clipping on supplementary alignments, those records carry a full copy of SEQ -- which is
# why 5,961 reads produced a 593 MB BAM upstream, and why the R step then loaded all of it into
# memory. Nothing downstream consumes them: processBamFn() indexes mapStats$primary throughout,
# and mapStats$secondary / $suppplemental feed only unused locals. Unmapped reads are NOT
# dropped (flag 0x4 is not in 0x900), so the read-length histogram is unaffected.
# Full flag accounting is preserved by teeing the unfiltered stream through flagstat.
# See docs/PATCHES.md.
minimap2 -t "$THREADS" -ax map-ont -Y "$REF" "$FASTQ" \
    | tee >(samtools flagstat - > "$FLAGSTAT") \
    | samtools view -u -F 0x900 - \
    | samtools sort -@ "$THREADS" -o "$OUT_BAM" -

samtools index "$OUT_BAM"

# The flagstat above runs in a process substitution, i.e. asynchronously, so it can still be
# writing when the pipeline returns. Give it a moment rather than racing it.
for _ in $(seq 1 100); do [[ -s "$FLAGSTAT" ]] && break; sleep 0.1; done

samtools quickcheck "$OUT_BAM"

echo "[align] done -> $OUT_BAM"
echo "[align] flagstat (before the secondary/supplementary filter):"
cat "$FLAGSTAT"
