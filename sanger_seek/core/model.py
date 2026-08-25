"""Core data model: project, samples, reads, alignments, variants."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .dna import complement, revcomp


@dataclass
class Config:
    """Analysis thresholds (kept explicit so they are auditable)."""

    trim_cutoff: float = 0.05      # Mott error-probability cutoff
    mixed_ratio: float = 0.30      # secondary/primary peak ratio -> mixed call
    high_qual: int = 30            # >= : strong support
    low_qual: int = 20             # <  : weak support
    near_end: int = 10             # bases from trim boundary considered fragile
    max_align_frac_dist: float = 0.40  # edit distance / read length above which alignment is rejected


@dataclass
class TraceData:
    """Raw chromatogram content of one AB1 file (never modified after load)."""

    channels: dict[str, np.ndarray]   # base -> analyzed trace (DATA9-12)
    calls: str                        # PBAS base calls
    quals: np.ndarray                 # PCON, uint8 per call
    ploc: np.ndarray                  # PLOC, trace sample index per call
    metadata: dict[str, str] = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return max((len(c) for c in self.channels.values()), default=0)


@dataclass
class PeakMetrics:
    """Per-called-base signal metrics (original AB1 call order)."""

    called: str
    primary_h: float
    secondary_h: float
    ratio: float                      # secondary / primary
    secondary_base: str
    noise: float


@dataclass
class DiscrepancyReport:
    """Imported .seq calls vs AB1-embedded calls. Both are retained."""

    count: int
    positions: list[int]              # PBAS coordinates (only when comparable)
    note: str = ""


@dataclass
class ReadAlignment:
    """Read aligned to the reference, in reference orientation."""

    orientation: str                  # 'F' or 'R'
    edit_distance: int
    identity: float
    ref_start: int                    # 0-based
    ref_end: int                      # exclusive
    read_start: int                   # oriented read coords, 0-based
    read_end: int                     # exclusive
    ops: list[tuple[str, int]]        # ('=', 'X', 'I', 'D') runs
    read_to_ref: np.ndarray           # oriented read index -> ref pos (-1 outside/inserted)
    ref_to_read: np.ndarray           # (refpos - ref_start) -> oriented read index (-1 deleted)

    def read_index_at(self, refpos: int) -> int | None:
        """Oriented read index aligned at refpos, else None."""
        if not (self.ref_start <= refpos < self.ref_end):
            return None
        i = int(self.ref_to_read[refpos - self.ref_start])
        return i if i >= 0 else None


@dataclass
class Read:
    id: str
    label: str
    ab1_path: Optional[str] = None
    seq_path: Optional[str] = None
    trace: Optional[TraceData] = None
    seq_imported: Optional[str] = None
    orientation_hint: Optional[str] = None

    # Derived during analysis
    calls: str = ""                   # authoritative calls (PBAS if trace, else .seq)
    quals: Optional[np.ndarray] = None
    trim: tuple[int, int] = (0, 0)    # [start, end) on original call order
    orientation: str = "?"            # 'F' | 'R' | '?'
    alignment: Optional[ReadAlignment] = None
    peaks: Optional[list[PeakMetrics]] = None
    discrepancies: Optional[DiscrepancyReport] = None
    error: str = ""

    @property
    def n(self) -> int:
        return len(self.calls)

    @property
    def is_reverse(self) -> bool:
        return self.orientation == "R"

    def orig_index(self, oriented_i: int) -> int:
        """Map an index in reference-oriented read coords to original AB1 coords."""
        return self.n - 1 - oriented_i if self.is_reverse else oriented_i

    def oriented_calls(self) -> str:
        return revcomp(self.calls) if self.is_reverse else self.calls

    def oriented_base(self, oriented_i: int) -> str:
        b = self.calls[self.orig_index(oriented_i)]
        return complement(b) if self.is_reverse else b

    def qual_at(self, oriented_i: int) -> Optional[int]:
        if self.quals is None:
            return None
        return int(self.quals[self.orig_index(oriented_i)])

    def peak_at(self, oriented_i: int) -> Optional[PeakMetrics]:
        if self.peaks is None:
            return None
        return self.peaks[self.orig_index(oriented_i)]

    def dist_from_trim_end(self, oriented_i: int) -> int:
        """Distance (bases) from the nearest trim boundary; <0 means outside trim."""
        oi = self.orig_index(oriented_i)
        s, e = self.trim
        return min(oi - s, e - 1 - oi)


@dataclass
class ReadEvidence:
    read_id: str
    read_label: str
    orientation: str
    covered: bool
    call: Optional[str]               # base(s) observed; '-' for deletion
    qual: Optional[int]
    supports: Optional[bool]          # None = ambiguous / no-call
    ratio: Optional[float]            # secondary/primary peak ratio at site
    secondary_base: Optional[str]
    oriented_index: Optional[int]
    dist_from_end: Optional[int]


@dataclass
class Variant:
    id: str
    kind: str                         # 'snv' | 'ins' | 'del'
    ref_pos: int                      # 0-based; ins = insertion before ref_pos
    ref_bases: str                    # '' for ins
    alt_bases: str                    # '' for del
    mixed: bool = False
    alt_fraction: Optional[float] = None
    evidence: list[ReadEvidence] = field(default_factory=list)
    strand_status: str = "one-read"   # 'both' | 'single' | 'discordant' | 'one-read'
    trace_quality: str = "clean"      # 'clean' | 'noisy'
    confidence: str = "review"        # 'high' | 'review' | 'low'
    cdna: Optional[str] = None
    protein: Optional[str] = None
    effect: Optional[str] = None      # missense/nonsense/synonymous/frameshift/...
    gene: Optional[str] = None
    control_status: str = "unavailable"  # present | absent | not-covered | unavailable

    @property
    def g_label(self) -> str:
        p = self.ref_pos + 1
        if self.kind == "snv":
            return f"g.{p}{self.ref_bases}>{self.alt_bases}"
        if self.kind == "del":
            if len(self.ref_bases) == 1:
                return f"g.{p}del{self.ref_bases}"
            return f"g.{p}_{p + len(self.ref_bases) - 1}del"
        return f"g.{p - 1}_{p}ins{self.alt_bases}"

    @property
    def label(self) -> str:
        return self.cdna or self.g_label


@dataclass
class CDSFeature:
    gene: str
    product: str
    strand: int                       # +1 / -1
    codon_start: int                  # 1..3
    parts: list[tuple[int, int]]      # genomic, 0-based half-open, file order
    genomic_order: np.ndarray         # transcription-order genomic indices
    cds_seq: str                      # spliced, stranded nucleotide sequence
    _g2c: Optional[dict[int, int]] = field(default=None, repr=False)

    def cds_index(self, gpos: int) -> Optional[int]:
        if self._g2c is None:
            self._g2c = {int(g): i for i, g in enumerate(self.genomic_order)}
        return self._g2c.get(gpos)

    def contains(self, gpos: int) -> bool:
        return any(s <= gpos < e for s, e in self.parts)


@dataclass
class Reference:
    name: str
    seq: str
    path: str = ""
    source: str = "fasta"             # 'fasta' | 'genbank' | 'seq'
    cds: list[CDSFeature] = field(default_factory=list)
    description: str = ""

    @property
    def n(self) -> int:
        return len(self.seq)

    def feature_at(self, gpos: int) -> Optional[CDSFeature]:
        for f in self.cds:
            if f.contains(gpos):
                return f
        return None


@dataclass
class Sample:
    key: str
    name: str
    reads: list[Read] = field(default_factory=list)
    variants: list[Variant] = field(default_factory=list)
    analyzed: bool = False
    error: str = ""

    def read_by_id(self, read_id: str) -> Optional[Read]:
        for r in self.reads:
            if r.id == read_id:
                return r
        return None

    @property
    def forward_reads(self) -> list[Read]:
        return [r for r in self.reads if r.orientation == "F"]

    @property
    def reverse_reads(self) -> list[Read]:
        return [r for r in self.reads if r.orientation == "R"]


@dataclass
class Project:
    reference: Optional[Reference] = None
    wt_control: Optional[Sample] = None
    samples: list[Sample] = field(default_factory=list)
    config: Config = field(default_factory=Config)
    path: str = ""

    def sample_by_key(self, key: str) -> Optional[Sample]:
        for s in self.samples:
            if s.key == key:
                return s
        return None
