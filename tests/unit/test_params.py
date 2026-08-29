import pytest
from pydantic import ValidationError

from telomers.params import TELOTAGS, TelomereRParams, TeloNPParams


def test_known_barcode_resolves_to_its_sequence():
    p = TelomereRParams(sample_label="s", barcode_name="NB65uq")
    assert p.barcode_seq == "TTCTCAGTCTTCCTCCAGACAAGG"


def test_unknown_barcode_names_the_alternatives():
    """Upstream failed with a bare "object 'NB99uq' not found" from inside R."""
    with pytest.raises(ValidationError) as exc:
        TelomereRParams(sample_label="s", barcode_name="NB99uq")
    message = str(exc.value)
    assert "NB65uq" in message and "barcode_seq" in message


def test_literal_barcode_bypasses_the_known_set():
    p = TelomereRParams(sample_label="s", barcode_name="custom", barcode_seq="acgtACGT")
    assert p.barcode_seq == "ACGTACGT"


def test_barcode_must_be_dna():
    with pytest.raises(ValidationError):
        TelomereRParams(sample_label="s", barcode_seq="ACGTX")


@pytest.mark.parametrize("bad", ["", "has space", "semi;colon", "a" * 65])
def test_sample_label_is_constrained(bad):
    """It ends up in filenames and a shell command line."""
    with pytest.raises(ValidationError):
        TeloNPParams(sample_label=bad)


def test_all_telotags_are_24bp_dna():
    assert len(TELOTAGS) == 18
    for name, seq in TELOTAGS.items():
        assert len(seq) == 24, name
        assert set(seq) <= set("ACGT"), name
