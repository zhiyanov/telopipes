"""Per-pipeline run parameters.

Hand-written Pydantic models rather than models generated from a spec file. With two pipelines
the generator would be more code than the thing it generates, and these carry domain detail --
the barcode enum-or-sequence rule, the note that chr_arm_ln comes from the reference -- that a
generic schema would lose. The web form is rendered from `model_json_schema()`, so the model
stays the single source of truth.
"""
from __future__ import annotations

import re
from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

#: The 24 bp TeloTag barcodes defined in telobp_functions.R, extracted from that file so the
#: two cannot drift. Upstream resolves these by name with get(); we pass the sequence itself
#: (see docs/PATCHES.md P5) but keep the names because they label the output columns.
TELOTAGS: dict[str, str] = {
    "NB01uq": "CACAAAGACACCGACAACTTTCTT",
    "NB02uq": "ACAGACGACTACAAACGGAATCGA",
    "NB10uq": "GAGAGGACAAAGGTTTCAACGCTT",
    "NB12uq": "TCCGATTCTGCTTCTTTCTACCTG",
    "NB13uq": "AGAACGACTTCCATACTCGTGTGA",
    "NB15uq": "AGGTCTACCTCGCTAACACCACTG",
    "NB16uq": "CGTCAACTGACAGTGGTTCGTACT",
    "NB19uq": "GTTCCTCGTGCAGTGTCAAGAGAT",
    "NB20uq": "TTGCGTCCTGTTACGAGAACTCAT",
    "NB50uq": "ATGGACTTTGGTAACTTCCTGCGT",
    "NB65uq": "TTCTCAGTCTTCCTCCAGACAAGG",
    "NB66uq": "CCGATCCTTGTGGCTTCTAACTTC",
    "NB67uq": "GTTTGTCATACTCGTGTGCTCACC",
    "NB68uq": "GAATCTAAGCAAACACGAAGGTGG",
    "NB69uq": "TACAGTCCGAGCCTCATGTGATCT",
    "NB70uq": "ACCGAGATCCTACGAATGGAGTGT",
    "NB72uq": "TAGCTGACTGTCTTCCATACCGAC",
    "NB88uq": "TCTTCTACTACCGATCCGAAGCAG",
}

_SAMPLE_RE = re.compile(r"^[A-Za-z0-9._-]{1,64}$")
_DNA_RE = re.compile(r"^[ACGT]+$")

SampleLabel = Annotated[
    str,
    Field(description="Label used in output filenames and plot titles",
          pattern=_SAMPLE_RE.pattern),
]


class TeloNPParams(BaseModel):
    """Pipeline B. Reference-free, so there is very little to configure."""

    sample_label: SampleLabel
    workers: int = Field(4, ge=1, le=32, description="Parallel scoring processes")
    chunk_reads: int = Field(
        2000, ge=100,
        description="Reads held in memory per dispatch round; bounds peak memory",
    )


class TelomereRParams(BaseModel):
    """Pipeline A. Mapping-based, so the reference constrains several of these."""

    sample_label: SampleLabel

    barcode_name: str = Field(
        "NB65uq",
        description="TeloTag name. Labels the output columns; also supplies the sequence "
                    "unless barcode_seq is given.",
    )
    barcode_seq: str | None = Field(
        None,
        description="Literal 24 bp TeloTag sequence, for barcodes outside the known set.",
    )

    #: NOT free to choose: the q-arm coordinate maths subtracts this, so a value that disagrees
    #: with the reference's per-record length silently corrupts q-arm lengths while leaving p
    #: arms perfect. The application derives it from the registered reference and shows it
    #: read-only. See docs/PLAN.md section 4.
    chr_arm_ln: int = Field(
        500_000, ge=1000,
        description="Per-record length of the cut reference. Read from the reference you "
                    "select, because a value that disagrees with it corrupts q-arm telomere "
                    "lengths while leaving p arms correct.",
    )

    read_len_min: int = Field(3000, ge=0, description="Drop reads shorter than this")
    aln_score_threshold: float = Field(
        20.0, ge=0, description="Minimum TeloTag alignment score"
    )
    map_start_threshold: int = Field(
        1000, ge=0,
        description="Maximum distance of the mapping start from the subtelomere boundary; "
                    "guards against interstitial telomere repeats",
    )
    tag_search_window: int = Field(
        300, ge=50, description="Bases at each read end searched for the barcode"
    )
    threads: int = Field(4, ge=1, le=32)

    @field_validator("barcode_seq")
    @classmethod
    def _check_seq(cls, v: str | None) -> str | None:
        if v is None:
            return None
        v = v.strip().upper()
        if not _DNA_RE.match(v):
            raise ValueError("barcode_seq must contain only A, C, G and T")
        return v

    @model_validator(mode="after")
    def _resolve_barcode(self):
        if self.barcode_seq is None:
            if self.barcode_name not in TELOTAGS:
                raise ValueError(
                    f"unknown barcode {self.barcode_name!r}. Known TeloTags: "
                    f"{', '.join(sorted(TELOTAGS))}. "
                    "Pass an arbitrary barcode with barcode_seq instead."
                )
            object.__setattr__(self, "barcode_seq", TELOTAGS[self.barcode_name])
        return self


PARAMS_MODELS: dict[str, type[BaseModel]] = {
    "telomere-r": TelomereRParams,
    "telonp": TeloNPParams,
}
