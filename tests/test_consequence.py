from sanger_seek.core.consequence import annotate_variants
from sanger_seek.core.model import Variant
from sanger_seek.core.reference import load_reference
from sanger_seek.devtools.demogen import (
    DEL_1BP,
    HET_SNV,
    INS_POS,
    NONSENSE_SNV,
    SYN_SNV,
)


def _v(kind, pos, ref, alt):
    return Variant(id="t", kind=kind, ref_pos=pos, ref_bases=ref, alt_bases=alt)


def test_missense(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.gb")
    v = _v("snv", HET_SNV, "C", "T")
    annotate_variants([v], ref)
    assert v.cdna == "c.944C>T"
    assert v.protein == "p.Thr315Ile"
    assert v.effect == "missense"
    assert v.gene == "SEEK1"


def test_synonymous(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.gb")
    v = _v("snv", SYN_SNV, "G", "A")
    annotate_variants([v], ref)
    assert v.cdna == "c.1173G>A"
    assert v.effect == "synonymous"
    assert v.protein == "p.Leu391="


def test_nonsense(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.gb")
    v = _v("snv", NONSENSE_SNV, "G", "A")
    annotate_variants([v], ref)
    assert v.cdna == "c.750G>A"
    assert v.effect == "nonsense"
    assert v.protein == "p.Trp250Ter"


def test_frameshift_deletion(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.gb")
    v = _v("del", DEL_1BP, "A", "")
    annotate_variants([v], ref)
    assert v.cdna == "c.1012delA"
    assert v.effect == "frameshift"
    assert v.protein == "p.Ser338fs"


def test_inframe_insertion(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.gb")
    v = _v("ins", INS_POS, "", "CTG")
    annotate_variants([v], ref)
    assert v.cdna == "c.600_601insCTG"
    assert v.effect == "in-frame insertion"


def test_noncoding(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.gb")
    v = _v("snv", 50, ref.seq[50], "A" if ref.seq[50] != "A" else "G")
    annotate_variants([v], ref)
    assert v.effect == "non-coding"
    assert v.cdna is None
    assert v.label.startswith("g.51")
