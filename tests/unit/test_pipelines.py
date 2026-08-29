from telomers import pipelines


def test_both_pipelines_load():
    reg = pipelines.registry()
    assert set(reg) == {"telomere-r", "telonp"}


def test_capabilities_drive_the_ui():
    """Pipeline B is reference-free and therefore cannot assign reads to chromosomes.
    The per-arm view must be absent rather than empty."""
    assert pipelines.get("telomere-r").has("per_arm")
    assert not pipelines.get("telonp").has("per_arm")


def test_reference_requirement():
    assert pipelines.get("telomere-r").requires_reference
    assert not pipelines.get("telonp").requires_reference


def test_every_pipeline_ships_a_snakefile_and_params_model():
    for pipeline in pipelines.registry().values():
        assert pipeline.snakefile.is_file()
        assert pipeline.params_model is not None
