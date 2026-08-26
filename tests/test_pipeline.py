"""End-to-end: demo dataset -> pairing -> analysis -> variants."""

import pytest

from sanger_seek.core.export import export_variants_csv
from sanger_seek.core.model import Config, Project
from sanger_seek.core.pairing import scan_paths
from sanger_seek.core.pipeline import analyze_sample, build_samples_from_scan
from sanger_seek.core.projectio import load_project, save_project
from sanger_seek.core.reference import load_reference
from sanger_seek.devtools.demogen import DEL_1BP, HET_SNV, INS_POS, NONSENSE_SNV, SYN_SNV


@pytest.fixture(scope="module")
def analyzed(demo_dir):
    scan = scan_paths([demo_dir])
    assert len(scan.references) == 2  # .gb and .fasta
    gb = next(p for p in scan.references if p.suffix == ".gb")
    reference = load_reference(gb)
    samples = build_samples_from_scan(scan)
    cfg = Config()
    for s in samples:
        analyze_sample(s, reference, cfg)
    return {s.key: s for s in samples}, reference


def test_pairing(analyzed):
    samples, _ = analyzed
    assert set(samples) == {"Sample001", "Sample002", "Sample003"}
    s1 = samples["Sample001"]
    assert len(s1.reads) == 2
    for r in s1.reads:
        assert r.ab1_path and r.seq_path


def test_orientation_and_alignment(analyzed):
    samples, _ = analyzed
    s1 = samples["Sample001"]
    by_label = {r.label: r for r in s1.reads}
    assert by_label["Sample001_F"].orientation == "F"
    assert by_label["Sample001_R"].orientation == "R"
    for r in s1.reads:
        assert r.alignment is not None
        assert r.alignment.identity > 0.98
        # junk ends were trimmed before alignment
        ts, te = r.trim
        assert ts >= 20 and te <= r.n - 20


def test_seq_discrepancies_are_retained(analyzed):
    samples, _ = analyzed
    s1 = samples["Sample001"]
    f = next(r for r in s1.reads if r.label == "Sample001_F")
    assert f.discrepancies is not None
    assert f.discrepancies.count == 2
    assert f.discrepancies.positions == [400, 401]
    # the AB1 calls were NOT overwritten by the .seq
    assert f.calls == f.trace.calls
    r = next(r for r in s1.reads if r.label == "Sample001_R")
    assert r.discrepancies.count == 0


def test_sample001_variants(analyzed):
    samples, _ = analyzed
    s1 = samples["Sample001"]
    by_pos = {(v.kind, v.ref_pos): v for v in s1.variants}

    het = by_pos[("snv", HET_SNV)]
    assert het.alt_bases == "T"
    assert het.mixed
    assert het.cdna == "c.944C>T"
    assert het.protein == "p.Thr315Ile"
    assert het.strand_status == "both"
    assert het.confidence == "high"
    assert het.alt_fraction is not None and 0.25 < het.alt_fraction < 0.55

    dele = by_pos[("del", DEL_1BP)]
    assert dele.ref_bases == "A"
    assert dele.cdna == "c.1012delA"
    assert dele.effect == "frameshift"
    assert dele.strand_status == "both"
    assert dele.confidence == "high"

    syn = by_pos[("snv", SYN_SNV)]
    assert syn.cdna == "c.1173G>A"
    assert syn.effect == "synonymous"
    # reverse read is a noisy no-call there -> single-strand support, review
    assert syn.strand_status == "single"
    assert syn.confidence == "review"

    # no spurious high-confidence calls beyond the three engineered ones
    extra = [v for v in s1.variants if (v.kind, v.ref_pos) not in by_pos or v.confidence == "high"]
    assert {(v.kind, v.ref_pos) for v in extra if v.confidence == "high"} == {
        ("snv", HET_SNV), ("del", DEL_1BP),
    }


def test_sample002_clean(analyzed):
    samples, _ = analyzed
    s2 = samples["Sample002"]
    assert [v for v in s2.variants if v.confidence != "low"] == []


def test_sample003_variants(analyzed):
    samples, _ = analyzed
    s3 = samples["Sample003"]
    by_pos = {(v.kind, v.ref_pos): v for v in s3.variants}

    non = by_pos[("snv", NONSENSE_SNV)]
    assert non.cdna == "c.750G>A"
    assert non.effect == "nonsense"
    assert non.protein == "p.Trp250Ter"
    assert non.strand_status == "one-read"

    ins = by_pos[("ins", INS_POS)]
    assert ins.alt_bases == "CTG"
    assert ins.cdna == "c.600_601insCTG"
    assert ins.effect == "in-frame insertion"


def test_export_and_project_roundtrip(analyzed, demo_dir, tmp_path):
    samples, reference = analyzed
    csv_path = tmp_path / "variants.csv"
    n = export_variants_csv(csv_path, list(samples.values()), reference)
    assert n >= 5
    text = csv_path.read_text()
    assert "c.944C>T" in text and "p.Thr315Ile" in text

    project = Project(
        reference=reference,
        wt_control=samples["Sample002"],
        samples=[samples["Sample001"], samples["Sample003"]],
    )
    project.samples[0].reads[0].orientation_override = "F"
    ppath = tmp_path / "demo.sanger-seek.json"
    save_project(project, ppath)
    loaded, warnings = load_project(ppath)
    assert warnings == []
    assert loaded.reference is not None and loaded.reference.name == reference.name
    assert loaded.wt_control is not None and loaded.wt_control.key == "Sample002"
    assert {s.key for s in loaded.samples} == {"Sample001", "Sample003"}
    assert len(loaded.sample_by_key("Sample001").reads) == 2
    assert loaded.sample_by_key("Sample001").reads[0].orientation_override == "F"
