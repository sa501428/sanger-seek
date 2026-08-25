import numpy as np

from sanger_seek.core import align as al
from sanger_seek.core.dna import revcomp
from sanger_seek.core.model import Config, Read

REF = "ACGTACGGATCCTTAGCCGGAATTCGCGATATCGGCCATGCACGTGTTGACCTGA"


def _read(calls: str, hint=None) -> Read:
    r = Read(id="t", label="t", orientation_hint=hint)
    r.calls = calls
    r.quals = np.full(len(calls), 50, dtype=np.uint8)
    r.trim = (0, len(calls))
    return r


def test_exact_match_forward():
    read = _read(REF[5:45])
    aln = al.orient_and_align(read, REF, Config())
    assert aln is not None
    assert aln.orientation == "F"
    assert aln.ref_start == 5 and aln.ref_end == 45
    assert aln.edit_distance == 0
    assert aln.identity == 1.0
    # coordinate maps are mutually consistent
    for i in range(40):
        assert aln.read_to_ref[i] == 5 + i
        assert aln.ref_to_read[i] == i


def test_reverse_orientation_detected():
    read = _read(revcomp(REF[5:45]), hint=None)
    aln = al.orient_and_align(read, REF, Config())
    assert aln is not None
    assert aln.orientation == "R"
    assert aln.ref_start == 5 and aln.ref_end == 45


def test_snv_and_maps():
    seq = list(REF[5:45])
    seq[10] = "A" if seq[10] != "A" else "G"
    read = _read("".join(seq))
    aln = al.orient_and_align(read, REF, Config())
    assert aln is not None
    assert aln.edit_distance == 1
    assert ("X", 1) in aln.ops
    assert aln.read_to_ref[10] == 15


def test_deletion_in_read():
    seq = REF[5:20] + REF[23:45]  # 3bp deletion
    read = _read(seq)
    aln = al.orient_and_align(read, REF, Config())
    assert aln is not None
    # co-optimal paths may split the run; total deleted bases is what matters
    assert sum(n for op, n in aln.ops if op == "D") == 3
    assert aln.edit_distance == 3
    assert (aln.ref_to_read == -1).sum() == 3


def test_insertion_in_read():
    seq = REF[5:25] + "TTT" + REF[25:45]
    read = _read(seq)
    aln = al.orient_and_align(read, REF, Config())
    assert aln is not None
    assert ("I", 3) in aln.ops
    ins_idx = 20
    assert aln.read_to_ref[ins_idx] == -1 or aln.read_to_ref[ins_idx + 2] == -1


def test_iupac_ambiguity_is_match():
    seq = list(REF[5:45])
    orig = seq[12]
    seq[12] = {"A": "R", "G": "R", "C": "Y", "T": "Y"}[orig]
    read = _read("".join(seq))
    aln = al.orient_and_align(read, REF, Config())
    assert aln is not None
    assert aln.edit_distance == 0


def test_reference_can_align_inside_longer_read_without_overhang_penalty():
    read = _read("TTTTAAAACC" + REF + "GGAATTTTCC", hint="F")
    aln = al.orient_and_align(read, REF, Config())
    assert aln is not None
    assert aln.orientation == "F"
    assert aln.edit_distance == 0
    assert aln.identity == 1.0
    assert (aln.ref_start, aln.ref_end) == (0, len(REF))
    assert (aln.read_start, aln.read_end) == (10, 10 + len(REF))


def test_python_fallback_matches_edlib():
    query = REF[5:20] + "T" + REF[20:40]
    raw_py = al._py_align(query, REF)
    assert raw_py is not None
    assert raw_py.edit_distance == 1
    if al.HAVE_EDLIB:
        raw_ed = al.align_infix(query, REF)
        assert raw_ed.edit_distance == raw_py.edit_distance
        assert raw_ed.ref_start == raw_py.ref_start


def test_reject_garbage():
    read = _read("GGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCAAAAATTTTTGGGGGCCCCC")
    aln = al.orient_and_align(read, "ATCGATCGATCG" * 20, Config())
    assert aln is None
