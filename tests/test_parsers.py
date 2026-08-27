import pytest

from sanger_seek.core.pairing import classify, scan_paths, split_direction
from sanger_seek.core.seqfile import parse_mutation_surveyor_text, parse_seq_text
from sanger_seek.core.phd import parse_phd_text


def test_parse_plain_seq():
    assert parse_seq_text("ACGT\nacgt\n") == "ACGTACGT"


def test_parse_fasta_headered_seq():
    assert parse_seq_text(">sample_F exported\nACG T\n120 ACGT\n") == "ACGTACGT"


def test_parse_seq_rejects_garbage():
    with pytest.raises(ValueError):
        parse_seq_text("hello world this is not dna")


def test_parse_mutation_surveyor_metadata_translation_and_numbered_dna():
    text = """\
/Gene = "BRCA2";
/Exon_And_Note = "coding reference";
/Reading Frame (1,2,3) = 2;
/CDS = 2..24;
/Amplicon Id = "";
/Translation = "    1 LLEICLKLVG CKMKKGLSSS
   21 ACGT ACGT
 1121 VKEISDIVQR XQ";

        1  ATGTTGGAGA TCTGCCTGAA GACCT
       26  NRYKMSWBDH V
"""
    parsed = parse_mutation_surveyor_text(text, name="BRCA2.seq")
    assert parsed.gene == "BRCA2"
    assert parsed.metadata["CDS"] == "2..24"
    assert parsed.metadata["Reading Frame (1,2,3)"] == "2"
    assert "1121 VKEISDIVQR XQ" in parsed.translation
    # The numbered translation line containing ACGT is metadata, not DNA.
    assert parsed.sequence == "ATGTTGGAGATCTGCCTGAAGACCTNRYKMSWBDHV"
    assert parse_seq_text(text) == parsed.sequence


def test_split_direction():
    assert split_direction("Patient001_F") == ("Patient001", "F")
    assert split_direction("Patient001-REV") == ("Patient001", "R")
    assert split_direction("s14 forward") == ("s14", "F")
    assert split_direction("Patient001") == ("Patient001", None)
    # trailing direction-ish letters that are part of the name stay put
    assert split_direction("BRAF")[1] is None


def test_classify():
    assert classify("a/b/x.ab1") == "ab1"
    assert classify("x.seq") == "seq"
    assert classify("x.fasta") == "fasta"
    assert classify("x.gbk") == "genbank"
    assert classify("x.phd") == "phd"
    assert classify("x.ab1.phd.1") == "phd"
    assert classify("x.pdf") is None


def test_scan_groups_reads(tmp_path):
    for name in [
        "P1_F.ab1", "P1_F.seq", "P1_R.ab1", "P1_R.seq",
        "P2_F.ab1", "P1_F.ab1.phd.1", "ref.gb", "notes.txt",
    ]:
        (tmp_path / name).write_bytes(b"x")
    result = scan_paths([tmp_path])
    assert set(result.samples) == {"P1", "P2"}
    assert len(result.samples["P1"]) == 2
    p1f = result.samples["P1"]["p1_f"]
    assert p1f.ab1 and p1f.seq and p1f.hint == "F"
    assert p1f.phd and p1f.phd.name == "P1_F.ab1.phd.1"
    assert [p.name for p in result.references] == ["ref.gb"]


def test_parse_phd_calls_qualities_and_locations():
    phd = parse_phd_text("""\
BEGIN_SEQUENCE read
BEGIN_DNA
a 42 10
C 18 21 extra
n 5 33
END_DNA
END_SEQUENCE
""")
    assert phd.calls == "ACN"
    assert phd.quals.tolist() == [42, 18, 5]
    assert phd.ploc.tolist() == [10, 21, 33]
