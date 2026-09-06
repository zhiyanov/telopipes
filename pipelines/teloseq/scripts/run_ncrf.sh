#!/usr/bin/env bash
# Telo-seq step 1: find telomeric repeat tracts in every read.
#
# Usage: run_ncrf.sh <reads.fastq[.gz]> <out_prefix>
set -euo pipefail
FASTQ="$1"
OUT="$2"

echo "[teloseq] FASTQ -> FASTA"
if [[ "$FASTQ" == *.gz ]]; then
    zcat "$FASTQ"
else
    cat "$FASTQ"
fi | awk 'NR%4==1{printf(">%s\n", substr($1,2))} NR%4==2{print}' > "${OUT}.fasta"

# Parameters are Telo-seq's, from its config.yml: motif telomeric:TTAGGG, minlength 50,
# minmratio 0.85, --positionalevents --scoring=nanopore --stats=events.
echo "[teloseq] NCRF repeat finding"
NCRF telomeric:TTAGGG \
        --minlength=50 --minmratio=0.85 \
        --positionalevents --scoring=nanopore --stats=events \
    < "${OUT}.fasta" > "${OUT}.ncrf"

echo "[teloseq] NCRF -> per-alignment summary"
python3 /opt/NCRF/ncrf_cat.py "${OUT}.ncrf" \
    | python3 /opt/NCRF/ncrf_summary.py > "${OUT}.ncrf.summary"

echo "[teloseq] $(( $(wc -l < "${OUT}.ncrf.summary") - 1 )) alignments -> ${OUT}.ncrf.summary"
