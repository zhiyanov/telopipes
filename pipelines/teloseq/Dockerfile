# Telo-seq analysis, reference-free (Schmidt et al., Nat Commun 2024).
#
# Identifies telomeric repeat tracts in each read with Noise Cancelling Repeat Finder, then
# derives a per-read telomere length from them. NCRF never sees a reference, so this variant
# reports a bulk distribution only; pipelines/teloseq-mapped adds chromosome assignment.
#
# Scope note: the published method basecalls raw pod5 signal with a custom Bonito model. We
# start from standard-basecalled FASTQ, so this is Telo-seq's *analysis* applied to ordinary
# reads, not the paper's pipeline end to end.
FROM docker.io/library/python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
    && rm -rf /var/lib/apt/lists/*

# NCRF v1.01.00 — the version Telo-seq pins. The C sources are portable (no x86 intrinsics),
# so this builds natively on arm64. -Werror is stripped because newer gcc warns on code that
# was clean in 2019.
RUN git clone --branch v1.01.00 --depth 1 \
        https://github.com/makovalab-psu/NoiseCancellingRepeatFinder.git /opt/NCRF \
 && sed -i 's/-Werror//g' /opt/NCRF/Makefile \
 && make -C /opt/NCRF \
 && test -x /opt/NCRF/NCRF

# NCRF's helper scripts are Python 2 and the base image has no python2 at all. Convert them
# rather than installing an end-of-life interpreter. Only the three we actually use are
# converted -- the rest of NCRF's toolbox has more Python 2 in it (xrange and friends) and
# converting code we never call would be unverifiable churn.
#
# lib2to3 fixes the print statements but not two things: the removed `file()` builtin, and
# string.maketrans moving onto str. Both are patched explicitly below.
RUN python3 -m lib2to3 -w -n --no-diffs /opt/NCRF/ncrf_cat.py /opt/NCRF/ncrf_summary.py \
        /opt/NCRF/ncrf_parse.py >/dev/null \
 && sed -i 's/= file(/= open(/g' /opt/NCRF/ncrf_cat.py /opt/NCRF/ncrf_summary.py \
 && sed -i -e 's/^from string import maketrans$/# python 3: maketrans lives on str/' \
           -e 's/\bmaketrans(/str.maketrans(/g' /opt/NCRF/ncrf_parse.py \
 && python3 -c "import ast,sys; [ast.parse(open(f).read()) for f in sys.argv[1:]]" \
        /opt/NCRF/ncrf_cat.py /opt/NCRF/ncrf_summary.py /opt/NCRF/ncrf_parse.py

# Prove the converted toolchain actually works, at build time rather than on someone's data:
# a synthetic read with a known 240 bp telomeric tract must survive NCRF -> cat -> summary.
RUN set -eu; \
    export PYTHONPATH=/opt/NCRF; \
    python3 -c "print('>t1'); print('TTAGGG'*40 + 'ACGTTGCATTACGCATGGATCCTA'*6)" > /tmp/t.fa; \
    /opt/NCRF/NCRF telomeric:TTAGGG --minlength=50 --minmratio=0.85 --positionalevents \
        --scoring=nanopore --stats=events < /tmp/t.fa > /tmp/t.ncrf; \
    python3 /opt/NCRF/ncrf_cat.py /tmp/t.ncrf | python3 /opt/NCRF/ncrf_summary.py > /tmp/t.sum; \
    cat /tmp/t.sum; \
    awk -F"\t" 'NR==2 && $3=="t1" && ($5-$4)==240 && $10==240 {ok=1} END{exit !ok}' /tmp/t.sum; \
    rm -f /tmp/t.fa /tmp/t.ncrf /tmp/t.sum; \
    echo "NCRF toolchain OK"

ENV PATH="/opt/NCRF:${PATH}" PYTHONPATH="/opt/NCRF" MPLBACKEND=Agg

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY scripts/ /opt/telomers/
COPY Snakefile /opt/telomers/Snakefile
WORKDIR /work
