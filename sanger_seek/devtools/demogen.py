"""Synthetic demo dataset generator.

Creates a small reference gene (FASTA + annotated GenBank) and matched
.ab1/.seq read pairs containing known variants, so the app and tests have
fully local, reproducible data:

  Sample001  het c.944C>T (p.Thr315Ile), hom c.1012delA (frameshift),
             hom c.1173G>A (synonymous; noisy/no-call on the reverse read)
  Sample002  clean bidirectional control
  Sample003  hom c.750G>A (p.Trp250Ter) + in-frame c.600_601insCTG
             (forward read only)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core.dna import complement
from .abif_writer import Tag, tag_char, tag_pstring, tag_short, write_abif

UTR5 = 200
CDS_LEN = 1200
UTR3 = 100

_SAFE_CODONS = [
    "GCT", "GCC", "TGC", "GAT", "GAA", "TTC", "GGA", "CAT", "ATT", "AAA",
    "CTT", "ATG", "AAC", "CCG", "CAA", "CGT", "TCT", "ACA", "GTG", "TGG",
]


def build_reference_seq(rng: np.random.Generator) -> str:
    codons = [
        _SAFE_CODONS[int(rng.integers(len(_SAFE_CODONS)))] for _ in range(CDS_LEN // 3)
    ]
    codons[0] = "ATG"
    codons[-1] = "TAA"
    codons[199] = "CAT"   # ends in T: keeps the CTG insertion left-aligned at c.600_601
    codons[249] = "TGG"   # Trp250 -> TGA nonsense target (c.750G>A)
    codons[314] = "ACC"   # Thr315 -> ATC missense target (c.944C>T)
    codons[336] = "CTC"   # context so the deleted A below is unambiguous
    codons[337] = "AGT"   # c.1012 = 'A' -> 1bp deletion target
    codons[390] = "CTG"   # Leu391 -> CTA synonymous target (c.1173G>A)
    cds = "".join(codons)
    assert len(cds) == CDS_LEN

    bases = "ACGT"
    utr5 = "".join(bases[int(i)] for i in rng.integers(4, size=UTR5))
    utr3 = "".join(bases[int(i)] for i in rng.integers(4, size=UTR3))
    return utr5 + cds + utr3


def write_reference_files(genome: str, out_dir: Path) -> None:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(
        Seq(genome),
        id="SEEK1_ref",
        name="SEEK1",
        description="Synthetic demo reference for sanger-seek",
    )
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"
    cds_seq = Seq(genome[UTR5 : UTR5 + CDS_LEN])
    rec.features.append(
        SeqFeature(
            FeatureLocation(UTR5, UTR5 + CDS_LEN, strand=1),
            type="CDS",
            qualifiers={
                "gene": ["SEEK1"],
                "product": ["seek demo protein 1"],
                "codon_start": ["1"],
                "translation": [str(cds_seq.translate())[:-1]],
            },
        )
    )
    SeqIO.write([rec], str(out_dir / "SEEK1.gb"), "genbank")
    with open(out_dir / "SEEK1.fasta", "w") as fh:
        fh.write(">SEEK1_ref Synthetic demo reference for sanger-seek\n")
        for i in range(0, len(genome), 70):
            fh.write(genome[i : i + 70] + "\n")


@dataclass
class SiteSpec:
    call: str                    # base call written to PBAS (IUPAC allowed)
    primary: str                 # dominant trace channel
    secondary: str | None = None
    frac2: float = 0.0           # secondary channel fraction of total signal


def _rc_spec(s: SiteSpec) -> SiteSpec:
    return SiteSpec(
        call=complement(s.call),
        primary=complement(s.primary),
        secondary=complement(s.secondary) if s.secondary else None,
        frac2=s.frac2,
    )


def build_read_entries(
    genome: str,
    region: tuple[int, int],
    per_pos: dict[int, SiteSpec] | None = None,
    deletions: set[int] | None = None,
    insertions: dict[int, str] | None = None,
    reverse: bool = False,
) -> tuple[list[SiteSpec], list[int | None]]:
    """Entries (call/trace spec per base) + genomic position per entry."""
    per_pos = per_pos or {}
    deletions = deletions or set()
    insertions = insertions or {}
    entries: list[SiteSpec] = []
    gpos: list[int | None] = []
    for g in range(*region):
        if g in insertions:
            for b in insertions[g]:
                entries.append(SiteSpec(b, b))
                gpos.append(None)
        if g in deletions:
            continue
        spec = per_pos.get(g, SiteSpec(genome[g], genome[g]))
        entries.append(spec)
        gpos.append(g)
    if reverse:
        entries = [_rc_spec(s) for s in reversed(entries)]
        gpos = list(reversed(gpos))
    return entries, gpos


def synth_trace(
    entries: list[SiteSpec],
    quals: np.ndarray,
    noisy: list[bool],
    rng: np.random.Generator,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    spacing, sigma, margin = 12, 3.0, 40
    n = len(entries)
    n_samples = margin * 2 + n * spacing
    x = np.arange(n_samples, dtype=np.float64)
    channels = {b: np.zeros(n_samples) for b in "ACGT"}
    ploc = np.zeros(n, dtype=np.int64)

    def add_peak(base: str, center: float, amp: float) -> None:
        lo = max(int(center - 5 * sigma), 0)
        hi = min(int(center + 5 * sigma) + 1, n_samples)
        channels[base][lo:hi] += amp * np.exp(
            -((x[lo:hi] - center) ** 2) / (2 * sigma * sigma)
        )

    for i, s in enumerate(entries):
        c = margin + i * spacing + spacing // 2
        ploc[i] = c
        amp = float(rng.uniform(950, 1500))
        if s.secondary and s.frac2 > 0:
            add_peak(s.primary, c, amp * (1 - s.frac2))
            add_peak(s.secondary, c, amp * s.frac2)
        else:
            add_peak(s.primary, c, amp)
        if noisy[i]:
            others = [b for b in "ACGT" if b != s.primary]
            for b in rng.choice(others, size=2, replace=False):
                add_peak(str(b), c + float(rng.normal(0, 2.0)), amp * float(rng.uniform(0.25, 0.45)))

    for b in "ACGT":
        channels[b] += np.abs(rng.normal(20, 8, n_samples))
        channels[b] = np.clip(channels[b], 0, 30000)
    return {b: channels[b].astype(np.int16) for b in "ACGT"}, ploc


def make_quals(
    n: int, junk_head: int, junk_tail: int, noisy: list[bool], rng: np.random.Generator
) -> np.ndarray:
    q = np.full(n, 50.0)
    ramp = 40
    for i in range(n):
        edge = min(i, n - 1 - i)
        if edge < ramp:
            q[i] = 15 + (50 - 15) * edge / ramp
    q += rng.normal(0, 2.0, n)
    q[:junk_head] = 8
    if junk_tail:
        q[-junk_tail:] = 8
    for i, flag in enumerate(noisy):
        if flag:
            q[i] = rng.uniform(14, 18)
    return np.clip(q, 2, 60).astype(np.uint8)


def write_read(
    out_dir: Path,
    name: str,
    genome: str,
    region: tuple[int, int],
    per_pos: dict[int, SiteSpec] | None = None,
    deletions: set[int] | None = None,
    insertions: dict[int, str] | None = None,
    reverse: bool = False,
    noisy_span: tuple[int, int] | None = None,
    seq_edits: list[int] | None = None,
    rng: np.random.Generator | None = None,
) -> None:
    rng = rng or np.random.default_rng(0)
    entries, gpos = build_read_entries(genome, region, per_pos, deletions, insertions, reverse)

    junk = 25
    bases = "ACGT"
    head = [SiteSpec(bases[int(i)], bases[int(i)]) for i in rng.integers(4, size=junk)]
    tail = [SiteSpec(bases[int(i)], bases[int(i)]) for i in rng.integers(4, size=junk)]
    entries = head + entries + tail
    gpos = [None] * junk + gpos + [None] * junk

    noisy = [
        g is not None and noisy_span is not None and noisy_span[0] <= g < noisy_span[1]
        for g in gpos
    ]
    quals = make_quals(len(entries), junk, junk, noisy, rng)
    channels, ploc = synth_trace(entries, quals, noisy, rng)
    calls = "".join(s.call for s in entries)

    order = "GATC"  # FWO_ channel order
    tags: list[Tag] = [
        tag_char("FWO_", 1, order),
        tag_pstring("SMPL", 1, name),
        tag_pstring("MCHN", 1, "SyntheSeq"),
        tag_pstring("MODL", 1, "3730"),
        tag_pstring("DySN", 1, "Z-BigDyeV3"),
    ]
    for i, b in enumerate(order):
        tags.append(tag_short("DATA", 9 + i, channels[b]))
    for num in (1, 2):
        tags.append(tag_char("PBAS", num, calls))
        tags.append(tag_short("PLOC", num, ploc))
        tags.append(tag_char("PCON", num, bytes(quals)))
    write_abif(out_dir / f"{name}.ab1", tags)

    exported = list(calls)
    rotate = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for idx in seq_edits or []:
        if 0 <= idx < len(exported):
            exported[idx] = rotate.get(exported[idx], "A")
    with open(out_dir / f"{name}.seq", "w") as fh:
        fh.write(f">{name} exported base calls\n")
        s = "".join(exported)
        for i in range(0, len(s), 60):
            fh.write(s[i : i + 60] + "\n")


# Genomic coordinates (0-based) of the demo variants
HET_SNV = UTR5 + 943      # c.944 C>T
DEL_1BP = UTR5 + 1011     # c.1012 delA
SYN_SNV = UTR5 + 1172     # c.1173 G>A
NONSENSE_SNV = UTR5 + 749  # c.750 G>A (TGG->TGA)
INS_POS = UTR5 + 600      # insertion before this position: c.600_601insCTG


def generate_demo(out_dir: str | Path, seed: int = 42) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    genome = build_reference_seq(rng)
    assert genome[HET_SNV] == "C" and genome[DEL_1BP] == "A" and genome[SYN_SNV] == "G"
    assert genome[NONSENSE_SNV] == "G"
    write_reference_files(genome, out)

    # Sample001: het C>T + 1bp del on both strands; syn SNV no-call on reverse
    write_read(
        out, "Sample001_F", genome, (140, 1460),
        per_pos={
            HET_SNV: SiteSpec("Y", "C", "T", 0.40),
            SYN_SNV: SiteSpec("A", "A"),
        },
        deletions={DEL_1BP},
        seq_edits=[400, 401],  # simulated manual edits in the .seq export
        rng=rng,
    )
    write_read(
        out, "Sample001_R", genome, (300, 1460),
        per_pos={
            HET_SNV: SiteSpec("T", "T", "C", 0.40),
            SYN_SNV: SiteSpec("N", "A", "G", 0.45),
        },
        deletions={DEL_1BP},
        reverse=True,
        noisy_span=(SYN_SNV - 15, SYN_SNV + 16),
        rng=rng,
    )

    # Sample002: clean control
    write_read(out, "Sample002_F", genome, (180, 1300), rng=rng)
    write_read(out, "Sample002_R", genome, (180, 1300), reverse=True, rng=rng)

    # Sample003: nonsense SNV + in-frame insertion (single read)
    write_read(
        out, "Sample003_F", genome, (400, 1300),
        per_pos={NONSENSE_SNV: SiteSpec("A", "A")},
        insertions={INS_POS: "CTG"},
        rng=rng,
    )
    return out
