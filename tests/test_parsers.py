import pytest

from sanger_seek.core.pairing import classify, scan_paths, split_direction
from sanger_seek.core.seqfile import parse_seq_text


def test_parse_plain_seq():
    assert parse_seq_text("ACGT\nacgt\n") == "ACGTACGT"


def test_parse_fasta_headered_seq():
    assert parse_seq_text(">sample_F exported\nACG T\n120 ACGT\n") == "ACGTACGT"


def test_parse_seq_rejects_garbage():
    with pytest.raises(ValueError):
        parse_seq_text("hello world this is not dna")


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
    assert classify("x.pdf") is None


def test_scan_groups_reads(tmp_path):
    for name in [
        "P1_F.ab1", "P1_F.seq", "P1_R.ab1", "P1_R.seq",
        "P2_F.ab1", "ref.gb", "notes.txt",
    ]:
        (tmp_path / name).write_bytes(b"x")
    result = scan_paths([tmp_path])
    assert set(result.samples) == {"P1", "P2"}
    assert len(result.samples["P1"]) == 2
    p1f = result.samples["P1"]["p1_f"]
    assert p1f.ab1 and p1f.seq and p1f.hint == "F"
    assert [p.name for p in result.references] == ["ref.gb"]
