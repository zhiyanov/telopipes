from pathlib import Path

from telomers import validate

REFS = Path(__file__).resolve().parents[1] / "data" / "refs"


def test_good_reference_passes_and_derives_its_arm_length():
    info = validate.validate_reference(REFS / "ok.fasta")
    assert info.ok
    assert info.n_records == 8
    # Derived from the file, never taken from the user: pipeline A subtracts this from q-arm
    # coordinates, so a wrong value corrupts q arms while leaving p arms perfect.
    assert info.arm_length_bp == 2000
    assert info.haplotypes == ["MATERNAL"]


def test_unpadded_contigs_are_rejected():
    """chr1p does not error in the pipeline -- the arm regex still matches -- it just sorts
    wrongly everywhere, because the R's own zero-padding fix-up is dead code."""
    info = validate.validate_reference(REFS / "unpadded.fasta")
    assert not info.ok
    assert any("zero-padded" in e for e in info.errors)


def test_ragged_records_are_rejected_with_the_offender_named():
    info = validate.validate_reference(REFS / "ragged.fasta")
    assert not info.ok
    assert info.arm_length_bp is None
    message = " ".join(info.errors)
    assert "chr02q_MATERNAL" in message and "q-arm" in message


def test_fastq_sniffing(tmp_path):
    good = tmp_path / "a.fastq"
    good.write_text("@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n")
    info = validate.sniff_fastq(good)
    assert info.ok and info.n_reads == 2 and info.total_bases == 6


def test_fastq_rejects_a_fasta(tmp_path):
    bad = tmp_path / "a.fasta"
    bad.write_text(">r1\nACGT\n")
    assert not validate.sniff_fastq(bad).ok


def test_fastq_reports_a_truncated_record(tmp_path):
    bad = tmp_path / "a.fastq"
    bad.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n")
    info = validate.sniff_fastq(bad)
    assert not info.ok and "truncated" in info.errors[0]


def test_barcode_guessed_from_the_naming_convention():
    assert validate.guess_barcode("HG002.1.NB65uq.fastq") == "NB65uq"
    assert validate.guess_barcode("sample.fastq") is None
