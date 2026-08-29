#!/usr/bin/env bash
# Phase 0 / smoke test: run pipeline A's two steps on the reference dataset and compare the
# result against the published figures in reference/results/.
#
# This is the check that validates the native-arm64 port: the image moved from
# bioconductor/bioconductor_docker (amd64, R 4.3 / Bioc 3.18) to rocker/r-ver (R 4.4 / Bioc
# 3.19) with different distro minimap2 and samtools. If the numbers below still come out, the
# port is sound. See docs/ENVIRONMENT.md.
#
# Usage: scripts/reproduce-reference.sh [reference_project_dir] [workdir]
set -euo pipefail

REF_PROJECT="${1:-/home/antzhi/bio/reference}"
WORK="${2:-$(pwd)/data/phase0-repro}"
IMAGE="${IMAGE:-telomers/pipeline-a:0.1.0}"
THREADS="${THREADS:-8}"

READS="$REF_PROJECT/data/HG002.1.NB65uq.fastq"
REFFA="$REF_PROJECT/data/HG002.v0.7.TeloBP.cut.MATERNAL.fasta"
EXPNAME="HG002.1.NB65uq.MATERNAL"

for f in "$READS" "$REFFA"; do
    [[ -f "$f" ]] || { echo "missing input: $f" >&2; exit 1; }
done

mkdir -p "$WORK"

run_in_image() {
    podman run --rm --platform linux/arm64 \
        -v "$WORK:/work:z" \
        -v "$READS:/inputs/reads.fastq:ro,z" \
        -v "$REFFA:/refs/reference.fasta:ro,z" \
        -w /work "$IMAGE" "$@"
}

echo "== step 1/2: align =="
time run_in_image bash /opt/telomers/align.sh \
    /refs/reference.fasta /inputs/reads.fastq /work/aligned.sort.bam "$THREADS"

echo
echo "== step 2/2: R analysis =="
# NOTE: the vendored run_analysis.R is still unpatched, so it takes positional arguments and
# writes reads.table_<exp>.prim.Samp.<bc>_alnScTh.20.txt. Patching it is phase 1 work; the
# baseline has to reproduce before it is changed.
time run_in_image Rscript /opt/telomers/run_analysis.R \
    /work/aligned.sort.bam /work/r_analysis "$EXPNAME" 500000 NB65uq

echo
echo "== comparison against reference/results =="
OURS="$WORK/r_analysis/summary.overall_${EXPNAME}.csv"
THEIRS="$REF_PROJECT/results/r_analysis/summary.overall_${EXPNAME}.csv"
echo "--- ours ---";      cat "$OURS"   2>/dev/null || echo "(missing: $OURS)"
echo "--- reference ---"; cat "$THEIRS" 2>/dev/null || echo "(missing: $THEIRS)"
echo "--- diff (empty means identical) ---"
diff "$THEIRS" "$OURS" && echo "IDENTICAL"

echo
echo "--- per-chromosome diff ---"
diff "$REF_PROJECT/results/r_analysis/summary.by_chr_${EXPNAME}.csv" \
     "$WORK/r_analysis/summary.by_chr_${EXPNAME}.csv" && echo "IDENTICAL"
