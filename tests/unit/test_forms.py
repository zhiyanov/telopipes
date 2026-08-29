from telomers import forms
from telomers.params import TelomereRParams, TeloNPParams


def test_fields_come_from_the_model():
    names = {f.name for f in forms.describe(TelomereRParams)}
    assert names == set(TelomereRParams.model_fields)


def test_barcode_renders_as_a_select_of_known_tags():
    field = next(f for f in forms.describe(TelomereRParams) if f.name == "barcode_name")
    assert field.kind == "select"
    assert len(field.choices) == 18 and "NB65uq" in field.choices


def test_chr_arm_ln_is_derived_from_the_reference():
    """It must not be typeable: a value disagreeing with the reference corrupts q-arm lengths
    while leaving p arms correct, which looks like biology rather than a bug."""
    field = next(f for f in forms.describe(TelomereRParams) if f.name == "chr_arm_ln")
    assert field.derived_from == "reference"
    assert field.help  # the form shows why it cannot be edited


def test_bounds_are_carried_through():
    field = next(f for f in forms.describe(TeloNPParams) if f.name == "workers")
    assert (field.minimum, field.maximum) == (1, 32)


def test_blank_form_values_fall_back_to_model_defaults():
    values = forms.coerce(TeloNPParams, {"sample_label": "s", "workers": "", "chunk_reads": " "})
    assert values == {"sample_label": "s"}
    assert TeloNPParams(**values).workers == 4


def test_unknown_form_keys_are_ignored():
    values = forms.coerce(TeloNPParams, {"sample_label": "s", "csrf": "x", "pipeline": "telonp"})
    assert set(values) == {"sample_label"}
