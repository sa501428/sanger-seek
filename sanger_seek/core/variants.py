"""Per-read candidate difference calling.

Walks the read/reference alignment and emits candidates for SNVs, indels and
mixed (heterozygous-looking) positions. Mixed positions are detected from
IUPAC ambiguity calls AND from secondary peak ratios — an IUPAC call aligned
to one of its bases counts as a match ('=') during alignment, so every
aligned column is inspected, not just mismatches.
"""

from __future__ import annotations

from dataclasses import dataclass

from .dna import complement, expand
from .model import Config, Read


@dataclass(frozen=True)
class CandidateKey:
    kind: str        # 'snv' | 'ins' | 'del'
    ref_pos: int
    ref_bases: str
    alt_bases: str


@dataclass
class ReadCandidate:
    key: CandidateKey
    call: str                    # observed call (IUPAC possible; ins bases; '' for del)
    mixed: bool
    oriented_index: int | None   # anchor position in oriented read coords
    qual: int | None
    ratio: float | None
    secondary_base: str | None   # in reference orientation
    alt_fraction: float | None   # secondary/(primary+secondary) when mixed


def site_metrics(read: Read, oriented_i: int) -> tuple[int | None, float | None, str | None, float | None]:
    """(qual, ratio, secondary_base_in_ref_orientation, alt_fraction) at a read position."""
    qual = read.qual_at(oriented_i)
    pk = read.peak_at(oriented_i)
    if pk is None:
        return qual, None, None, None
    sec = complement(pk.secondary_base) if read.is_reverse else pk.secondary_base
    denom = pk.primary_h + pk.secondary_h
    frac = pk.secondary_h / denom if denom > 0 else None
    return qual, pk.ratio, sec, frac


def normalize_indel(kind: str, pos: int, bases: str, ref_seq: str) -> tuple[int, str]:
    """Left-align an indel so equivalent representations from different reads
    (or strands) share one canonical position and allele."""
    if kind == "del":
        k = len(bases)
        while pos > 0 and ref_seq[pos - 1] == ref_seq[pos + k - 1]:
            pos -= 1
        return pos, ref_seq[pos : pos + k]
    s = bases
    while pos > 0 and s and ref_seq[pos - 1] == s[-1]:
        s = ref_seq[pos - 1] + s[:-1]
        pos -= 1
    return pos, s


def read_candidates(read: Read, ref_seq: str, cfg: Config) -> list[ReadCandidate]:
    aln = read.alignment
    if aln is None:
        return []
    out: list[ReadCandidate] = []

    rc = aln.read_start
    fc = aln.ref_start
    ops = aln.ops
    for oi, (op, k) in enumerate(ops):
        if op in ("=", "X"):
            for t in range(k):
                cand = _snv_candidate(read, ref_seq, cfg, rc + t, fc + t)
                if cand is not None:
                    out.append(cand)
            rc += k
            fc += k
        elif op == "I":
            interior = 0 < oi < len(ops) - 1
            if interior:
                bases = "".join(read.oriented_base(rc + t) for t in range(k))
                quals = [read.qual_at(rc + t) for t in range(k)]
                quals = [q for q in quals if q is not None]
                _, ratio, sec, frac = site_metrics(read, rc)
                ins_pos, ins_bases = normalize_indel("ins", fc, bases, ref_seq)
                out.append(
                    ReadCandidate(
                        key=CandidateKey("ins", ins_pos, "", ins_bases),
                        call=bases,
                        mixed=False,
                        oriented_index=rc,
                        qual=min(quals) if quals else None,
                        ratio=ratio,
                        secondary_base=sec,
                        alt_fraction=frac,
                    )
                )
            rc += k
        elif op == "D":
            interior = 0 < oi < len(ops) - 1
            if interior:
                flank_quals = [
                    q
                    for q in (
                        read.qual_at(rc - 1) if rc - 1 >= 0 else None,
                        read.qual_at(rc) if rc < read.n else None,
                    )
                    if q is not None
                ]
                _, ratio, sec, frac = site_metrics(read, max(rc - 1, 0))
                del_pos, del_bases = normalize_indel("del", fc, ref_seq[fc : fc + k], ref_seq)
                out.append(
                    ReadCandidate(
                        key=CandidateKey("del", del_pos, del_bases, ""),
                        call="",
                        mixed=False,
                        oriented_index=max(rc - 1, 0),
                        qual=min(flank_quals) if flank_quals else None,
                        ratio=ratio,
                        secondary_base=sec,
                        alt_fraction=frac,
                    )
                )
            fc += k
    return out


def _snv_candidate(
    read: Read, ref_seq: str, cfg: Config, oriented_i: int, refpos: int
) -> ReadCandidate | None:
    ref_base = ref_seq[refpos]
    call = read.oriented_base(oriented_i)
    if call == "N":
        return None
    call_set = set(expand(call))
    if not call_set:
        return None
    qual, ratio, sec, frac = site_metrics(read, oriented_i)

    alts = call_set - {ref_base}
    if len(call_set) > 1:
        if not alts:
            return None
        # Ambiguity call: pick the alt supported by the secondary peak when possible.
        alt = sec if sec in alts else sorted(alts)[0]
        mixed = ref_base in call_set or (ratio is not None and ratio >= cfg.mixed_ratio)
        return ReadCandidate(
            key=CandidateKey("snv", refpos, ref_base, alt),
            call=call, mixed=mixed, oriented_index=oriented_i,
            qual=qual, ratio=ratio, secondary_base=sec, alt_fraction=frac,
        )
    if alts:
        # Plain mismatch; het if the reference base persists as a strong secondary.
        mixed = ratio is not None and ratio >= cfg.mixed_ratio and sec == ref_base
        return ReadCandidate(
            key=CandidateKey("snv", refpos, ref_base, call),
            call=call, mixed=mixed, oriented_index=oriented_i,
            qual=qual, ratio=ratio, secondary_base=sec, alt_fraction=frac,
        )
    # Call matches reference — check for an undercalled mixed peak.
    if (
        ratio is not None
        and sec is not None
        and ratio >= cfg.mixed_ratio
        and sec != ref_base
        and sec in "ACGT"
    ):
        return ReadCandidate(
            key=CandidateKey("snv", refpos, ref_base, sec),
            call=call, mixed=True, oriented_index=oriented_i,
            qual=qual, ratio=ratio, secondary_base=sec, alt_fraction=frac,
        )
    return None
