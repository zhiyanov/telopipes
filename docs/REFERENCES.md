# Choosing a reference

Two of the four pipelines need one; two never look at it.

| Pipeline | Reference | Reports |
|---|---|---|
| `telomere-r` | required, cut | per chromosome arm |
| `teloseq-mapped` | required, cut | per chromosome arm |
| `telonp` | none | bulk only |
| `teloseq` | none | bulk only |

## Why it has to be cut

`telomere-r` measures a telomere as the bases minimap2 had to soft-clip. That only works when
the reference has nothing to align the telomere against. Hand it a telomere-to-telomere
assembly — which genuinely carries `(TTAGGG)n` at every chromosome end — and the aligner
extends into the assembly's own telomere, so the clip becomes
*sample telomere − assembly telomere*, and any read shorter than the assembly's registers as
nothing at all.

So register the whole genome and let the application derive the cut form:

```bash
telomers references add --cut ~/Downloads/hg002v1.2.fasta.gz            # diploid, 92 arms
telomers references add --cut --haplotype MATERNAL hg002v1.2.fasta.gz   # haploid, 46 arms
```

Cutting HG002 v1.2 takes about two minutes and produces records named `chr01p_MATERNAL`,
`chr01q_MATERNAL`, … each exactly 500 kb. The naming matters and is normalised during cutting:
`posEndPQchrFn()` in the R pipeline zero-pads `chr1p` to `chr01p` on one line and overwrites the
result on the next, so an unpadded reference does not error — it just sorts wrongly throughout.

## Assembly version

Moving from HG002 v0.7 to v1.2 changes almost nothing in bulk — median 4,593.5 → 4,570.5 bp on
the same reads — and 36 of 46 arms agree within 50 bp. What it does fix is the acrocentrics:
`chr15p` loses 166 reads and `chr22p` gains 165, i.e. v1.2 separates two near-identical short
arms that v0.7 conflated, moving chr22p's median by +1,258 bp. `chr14p`, the shortest-telomere
result, is unchanged.

**Use the newest assembly.** The bulk answer is insensitive to it, and the per-arm answer is
better.

## Haploid or diploid

Measured on HG002.1.NB65uq with `telomere-r`, same reads, same assembly:

| | maternal (46 arms) | diploid (92 arms) |
|---|---|---|
| reads passing | 4,358 | 4,368 |
| median telomere | 4,593.5 bp | 4,570.5 bp |
| mapq ≥ 40 | **75.8%** | **21.1%** |
| mapq < 10 | 16.6% | 52.2% |
| mapq 0 | 0.7% | 9.0% |

**Bulk length does not care.** 23 bp apart, with the same number of reads passing.

**Allele assignment mostly does not survive.** With both haplotypes present, maternal and
paternal chromosome ends are near-identical, mapping confidence collapses, and reads divide
50.2% / 49.8% between the haplotypes — the signature of an aligner choosing between two equally
good placements, not of an allele call. Raising the mapq cutoff never makes the split decisive
(48/52 at ≥10, 46/54 at ≥40) and costs most of the reads.

Restricted to mapq ≥ 40, only **9 of 44 autosomal+X arms** retain ≥8 confident reads on *both*
alleles. Those nine do show large differences — `chr19q` 1,921 vs 5,510 bp, `chr04q` 4,160 vs
7,638 bp — which is interesting and is not something a haploid reference can show at all. But at
this depth it is a lead to follow, not a result: Karimian et al. reach allele-specific lengths
with a two-step diploid→haploid protocol (mapq 10, then 30), which this application does not
implement.

**Recommendation.** Use a haploid cut reference for per-arm telomere lengths — the numbers are
the same and the mapping is far more confident. Use the diploid one when allele-specific
differences are the actual question, and read the per-arm table together with the mapq column
rather than on its own.

## All four pipelines, diploid reference, HG002.1.NB65uq

| pipeline | passed / total | median | mean | arm×haplotype groups |
|---|---|---|---|---|
| `telomere-r` | 4,368 / 4,737 | 4,570 bp | 4,825 bp | 91 |
| `telonp` | 5,049 / 5,961 | 4,378 bp | 4,471 bp | — |
| `teloseq` | 2,521 / 4,234 | 4,665 bp | 5,774 bp | — |
| `teloseq-mapped` | 2,371 / 4,234 | 4,704 bp | 5,890 bp | 91 |

All four bulk medians fall within ~330 bp of each other, which is the reassuring part.

The two mapped pipelines agree well per arm: of 83 arm×haplotype groups with ≥10 reads in both,
the median difference is **110 bp** and 72 are within 500 bp. The disagreement is concentrated,
not diffuse — `chr19q` differs by **+21,088 bp on the maternal allele and +20,776 bp on the
paternal one**. That is NCRF's `end.max − start.min` spanning an interstitial telomeric repeat,
and seeing it on both alleles independently is what confirms the cause is the chromosome, not
the mapping. Reading `teloseq-mapped` output, check the `segments` column before trusting a
long tract.
