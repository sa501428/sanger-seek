from sanger_seek.core.dna import revcomp
from sanger_seek.core.reference import load_reference


def test_load_seq_as_explicit_reference(tmp_path):
    path = tmp_path / "wild_type.seq"
    path.write_text("ACGT ACGT\nNN\n")
    ref = load_reference(path)
    assert ref.name == "wild_type"
    assert ref.seq == "ACGTACGTNN"
    assert ref.source == "seq"
    assert ref.cds == []
from sanger_seek.devtools.demogen import CDS_LEN, UTR5


def test_load_genbank_reference(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.gb")
    assert ref.source == "genbank"
    assert len(ref.seq) == 1500
    assert len(ref.cds) == 1
    cds = ref.cds[0]
    assert cds.gene == "SEEK1"
    assert cds.strand == 1
    assert cds.parts == [(UTR5, UTR5 + CDS_LEN)]
    assert cds.cds_seq == ref.seq[UTR5 : UTR5 + CDS_LEN]
    assert cds.cds_seq.startswith("ATG") and cds.cds_seq.endswith("TAA")
    # genomic <-> cds mapping
    assert cds.cds_index(UTR5) == 0
    assert cds.cds_index(UTR5 + 943) == 943
    assert cds.cds_index(UTR5 - 1) is None
    assert ref.feature_at(UTR5 + 5) is cds
    assert ref.feature_at(10) is None


def test_load_genbank_record_with_seq_suffix(demo_dir, tmp_path):
    path = tmp_path / "exported-reference.seq"
    path.write_bytes((demo_dir / "SEEK1.gb").read_bytes())
    ref = load_reference(path)
    assert ref.source == "genbank"
    assert len(ref.cds) == 1


def test_load_fasta_reference(demo_dir):
    ref = load_reference(demo_dir / "SEEK1.fasta")
    assert ref.source == "fasta"
    assert ref.cds == []
    gb = load_reference(demo_dir / "SEEK1.gb")
    assert ref.seq == gb.seq


def test_minus_strand_cds(tmp_path):
    # Reference whose CDS is on the minus strand: ATG..TAA reverse-complemented
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    coding = "ATGGCTACCGAATAA"  # M A T E *
    genome = "TTTTT" + revcomp(coding) + "GGGGG"
    rec = SeqRecord(Seq(genome), id="mini", name="mini", description="")
    rec.annotations["molecule_type"] = "DNA"
    rec.features.append(
        SeqFeature(
            FeatureLocation(5, 5 + len(coding), strand=-1),
            type="CDS",
            qualifiers={"gene": ["MINI"]},
        )
    )
    path = tmp_path / "mini.gb"
    SeqIO.write([rec], str(path), "genbank")

    ref = load_reference(path)
    cds = ref.cds[0]
    assert cds.strand == -1
    assert cds.cds_seq == coding
    # last genomic base of the feature is the first cds base
    assert cds.cds_index(5 + len(coding) - 1) == 0
    assert cds.cds_index(5) == len(coding) - 1
