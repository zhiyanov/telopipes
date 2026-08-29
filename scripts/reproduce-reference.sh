#!/usr/bin/env bash
# Reproduce the reference results with the containerised pipelines.
#
# This is the check that validates the native-arm64 port. The reference was produced on amd64
# with bioconductor/bioconductor_docker (R 4.3 / Bioc 3.18); we run on arm64 with rocker/r-ver
# (R 4.4 / Bioc 3.19) and distro minimap2/samtools. Every output below is expected to be
# byte-identical anyway. See docs/ENVIRONMENT.md.
#
# Usage: scripts/reproduce-reference.sh [reference_project_dir] [workdir]
set -uo pipefail

REF_PROJECT="${1:-/home/antzhi/bio/reference}"
WORK="${2:-$(pwd)/data/reproduce}"
IMAGE_A="${IMAGE_A:-telomers/pipeline-a:0.1.0}"
IMAGE_B="${IMAGE_B:-telomers/pipeline-b:0.1.0}"
THREADS="${THREADS:-8}"

READS="$REF_PROJECT/data/HG002.1.NB65uq.fastq"
REFFA="$REF_PROJECT/data/HG002.v0.7.TeloBP.cut.MATERNAL.fasta"
SAMPLE="HG002.1.NB65uq"
EXPNAME="HG002.1.NB65uq.MATERNAL"
RES="$REF_PROJECT/results"

for f in "$READS" "$REFFA"; do
    [[ -f "$f" ]] || { echo "missing input: $f" >&2; exit 1; }
done

FAILURES=0
check() {  # check <label> <reference_file> <our_file>
    printf '  %-46s ' "$1"
    if [[ ! -f "$2" ]]; then echo "SKIP (no reference file)"; return; fi
    if [[ ! -f "$3" ]]; then echo "FAIL (not produced)"; FAILURES=$((FAILURES+1)); return; fi
    if diff -q "$2" "$3" >/dev/null 2>&1; then
        echo "identical"
    else
        echo "DIFFERS"; FAILURES=$((FAILURES+1)); diff "$2" "$3" | head -5 | sed 's/^/      /'
    fi
}

mkdir -p "$WORK/a" "$WORK/b"

# ---------------------------------------------------------------- pipeline A
# NOTE: run_analysis.R is still the unpatched upstream script, so it takes positional
# arguments and writes reads.table_<exp>.prim.Samp.<bc>_alnScTh.20.txt. Patching it is phase 1;
# the baseline has to reproduce before it is changed.
in_a() {
    podman run --rm --platform linux/arm64 \
        -v "$WORK/a:/work:z" \
        -v "$READS:/inputs/reads.fastq:ro,z" \
        -v "$REFFA:/refs/reference.fasta:ro,z" \
        -w /work "$IMAGE_A" "$@"
}

echo "== pipeline A, step 1/2: align =="
time in_a bash /opt/telomers/align.sh \
    /refs/reference.fasta /inputs/reads.fastq /work/aligned.sort.bam "$THREADS" 2>&1 | tail -14

echo
echo "== pipeline A, step 2/2: R analysis =="
time in_a Rscript /opt/telomers/run_analysis.R \
    /work/aligned.sort.bam /work/r_analysis "$EXPNAME" 500000 NB65uq 2>&1 | tail -6

# ---------------------------------------------------------------- pipeline B
# --shm-size is required: pandarallel moves DataFrame chunks through /dev/shm, which podman
# caps at 64 MB, and the default fails with ENOSPC on a machine with plenty of free disk.
echo
echo "== pipeline B: TeloNP =="
time podman run --rm --platform linux/arm64 --shm-size=2g \
    -v "$WORK/b:/work:z" \
    -v "$READS:/inputs/reads.fastq:ro,z" \
    -w /work "$IMAGE_B" \
    python /opt/telomers/run_telonp.py /inputs/reads.fastq /work --sample "$SAMPLE" \
    2>&1 | grep -vE '^\s*[0-9.]+%' | tail -20

# ---------------------------------------------------------------- comparison
echo
echo "=============================================================="
echo "  comparison against $RES"
echo "=============================================================="
echo "pipeline A:"
REFBAM="$RES/alignment/${EXPNAME}.sort.bam"
if [[ -f "$REFBAM" ]]; then
    printf '  %-46s ' "alignment records (headers excluded)"
    OURS=$(in_a sh -c 'samtools view /work/aligned.sort.bam | md5sum' 2>/dev/null | awk '{print $1}')
    THEIRS=$(podman run --rm --platform linux/arm64 -v "$RES/alignment:/ref:ro,z" "$IMAGE_A" \
             sh -c "samtools view /ref/$(basename "$REFBAM") | md5sum" 2>/dev/null | awk '{print $1}')
    if [[ -n "$OURS" && "$OURS" == "$THEIRS" ]]; then echo "identical  ($OURS)"
    else echo "DIFFERS"; echo "      ours=$OURS theirs=$THEIRS"; FAILURES=$((FAILURES+1)); fi
fi
check "reads.table (per read)" \
      "$RES/r_analysis/reads.table_${EXPNAME}.prim.Samp.NB65uq_alnScTh.20.txt" \
      "$WORK/a/r_analysis/reads.table_${EXPNAME}.prim.Samp.NB65uq_alnScTh.20.txt"
check "summary.overall" "$RES/r_analysis/summary.overall_${EXPNAME}.csv" \
      "$WORK/a/r_analysis/summary.overall_${EXPNAME}.csv"
check "summary.by_chr"  "$RES/r_analysis/summary.by_chr_${EXPNAME}.csv" \
      "$WORK/a/r_analysis/summary.by_chr_${EXPNAME}.csv"

echo "pipeline B:"
for f in "${SAMPLE}.telonp.all.csv" "${SAMPLE}.telonp.csv" "${SAMPLE}.telonp.summary.txt"; do
    check "$f" "$RES/telonp/$f" "$WORK/b/$f"
done

echo
if [[ "$FAILURES" -eq 0 ]]; then
    echo "RESULT: every output reproduces the reference byte-for-byte."
else
    echo "RESULT: $FAILURES check(s) did not match."
fi
exit "$FAILURES"
