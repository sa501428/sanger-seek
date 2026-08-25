"""Forward/reverse reconciliation: merge per-read candidates into sample
variants with strand status, trace quality and confidence.

For every candidate site each read contributes explicit evidence (covered?
call? quality? peak ratio?) so the UI can show exactly why a variant was
flagged and let the user audit it against the raw traces.
"""

from __future__ import annotations

from .dna import expand
from .model import Config, Read, ReadEvidence, Sample, Variant
from .variants import CandidateKey, ReadCandidate, read_candidates, site_metrics


def call_sample_variants(sample: Sample, ref_seq: str, cfg: Config) -> list[Variant]:
    per_read: dict[str, dict[CandidateKey, ReadCandidate]] = {}
    keys: dict[CandidateKey, None] = {}
    for read in sample.reads:
        cands = read_candidates(read, ref_seq, cfg)
        per_read[read.id] = {c.key: c for c in cands}
        for c in cands:
            keys.setdefault(c.key, None)

    variants: list[Variant] = []
    for key in keys:
        evidence = [
            _evidence_for(read, key, per_read.get(read.id, {}), ref_seq, cfg)
            for read in sample.reads
        ]
        if not any(e.supports for e in evidence):
            continue
        variants.append(_build_variant(key, evidence, per_read, sample, cfg))

    variants.sort(key=lambda v: (v.ref_pos, v.kind, v.alt_bases))
    return variants


def compare_sample_to_control(sample: Sample, control: Sample | None) -> None:
    """Label sample variants according to an independently aligned WT control."""
    if control is None or not control.analyzed:
        for variant in sample.variants:
            variant.control_status = "unavailable"
        return
    control_keys = {
        (v.kind, v.ref_pos, v.ref_bases, v.alt_bases) for v in control.variants
    }
    for variant in sample.variants:
        key = (variant.kind, variant.ref_pos, variant.ref_bases, variant.alt_bases)
        if key in control_keys:
            variant.control_status = "present"
            continue
        covered = False
        for read in control.reads:
            aln = read.alignment
            if aln is None:
                continue
            if variant.kind == "ins":
                covered = aln.ref_start <= variant.ref_pos - 1 and variant.ref_pos < aln.ref_end
            else:
                covered = aln.ref_start <= variant.ref_pos < aln.ref_end
            if covered:
                break
        variant.control_status = "absent" if covered else "not-covered"


def _evidence_for(
    read: Read,
    key: CandidateKey,
    cands: dict[CandidateKey, ReadCandidate],
    ref_seq: str,
    cfg: Config,
) -> ReadEvidence:
    aln = read.alignment
    ev = ReadEvidence(
        read_id=read.id,
        read_label=read.label,
        orientation=read.orientation,
        covered=False,
        call=None,
        qual=None,
        supports=None,
        ratio=None,
        secondary_base=None,
        oriented_index=None,
        dist_from_end=None,
    )
    if aln is None:
        return ev

    cand = cands.get(key)
    if key.kind == "snv":
        idx = aln.read_index_at(key.ref_pos)
        if not (aln.ref_start <= key.ref_pos < aln.ref_end):
            return ev
        ev.covered = True
        if idx is None:  # deleted in this read
            ev.call = "-"
            ev.supports = False
            return ev
        ev.oriented_index = idx
        ev.dist_from_end = read.dist_from_trim_end(idx)
        call = read.oriented_base(idx)
        ev.call = call
        ev.qual, ev.ratio, ev.secondary_base, _ = site_metrics(read, idx)
        if cand is not None:
            ev.supports = True
        elif call == "N":
            ev.supports = None
        elif key.alt_bases in expand(call):
            ev.supports = True
        else:
            ev.supports = False
        return ev

    # Indels: read must align across the affected region.
    if key.kind == "del":
        span_ok = aln.ref_start <= key.ref_pos - 1 and key.ref_pos + len(key.ref_bases) < aln.ref_end
    else:  # ins: junction between ref_pos-1 and ref_pos
        span_ok = aln.ref_start <= key.ref_pos - 1 and key.ref_pos < aln.ref_end
    if not span_ok:
        return ev
    ev.covered = True
    if cand is not None:
        ev.supports = True
        ev.call = "-" * len(key.ref_bases) if key.kind == "del" else cand.call
        ev.qual = cand.qual
        ev.ratio = cand.ratio
        ev.secondary_base = cand.secondary_base
        ev.oriented_index = cand.oriented_index
        if cand.oriented_index is not None:
            ev.dist_from_end = read.dist_from_trim_end(cand.oriented_index)
        return ev
    # No matching indel: report what the read shows across that region.
    ev.supports = False
    if key.kind == "del":
        calls = []
        quals = []
        for p in range(key.ref_pos, key.ref_pos + len(key.ref_bases)):
            idx = aln.read_index_at(p)
            if idx is None:
                calls.append("-")
            else:
                calls.append(read.oriented_base(idx))
                q = read.qual_at(idx)
                if q is not None:
                    quals.append(q)
        ev.call = "".join(calls)
        ev.qual = min(quals) if quals else None
        anchor = aln.read_index_at(key.ref_pos)
    else:
        ev.call = ""
        anchor = aln.read_index_at(key.ref_pos - 1)
        if anchor is not None:
            q = read.qual_at(anchor)
            ev.qual = q
    if anchor is not None:
        ev.oriented_index = anchor
        ev.dist_from_end = read.dist_from_trim_end(anchor)
        _, ev.ratio, ev.secondary_base, _ = site_metrics(read, anchor)
    return ev


def _build_variant(
    key: CandidateKey,
    evidence: list[ReadEvidence],
    per_read: dict[str, dict[CandidateKey, ReadCandidate]],
    sample: Sample,
    cfg: Config,
) -> Variant:
    covering = [e for e in evidence if e.covered]
    supporting = [e for e in covering if e.supports]
    contradicting = [e for e in covering if e.supports is False]

    sup_orients = {e.orientation for e in supporting if e.orientation in ("F", "R")}
    if len(covering) <= 1:
        strand_status = "one-read"
    elif supporting and contradicting:
        strand_status = "discordant"
    elif {"F", "R"} <= sup_orients:
        strand_status = "both"
    else:
        strand_status = "single"

    cands = [
        per_read[e.read_id][key]
        for e in supporting
        if key in per_read.get(e.read_id, {})
    ]
    mixed = any(c.mixed for c in cands)
    fracs = [c.alt_fraction for c in cands if c.mixed and c.alt_fraction is not None]
    alt_fraction = sum(fracs) / len(fracs) if fracs else None

    noisy = False
    for e in supporting:
        if e.qual is not None and e.qual < cfg.low_qual:
            noisy = True
        if e.ratio is not None and e.ratio > (0.90 if mixed else 0.50):
            noisy = True
        if e.dist_from_end is not None and e.dist_from_end < cfg.near_end:
            noisy = True
    trace_quality = "noisy" if noisy else "clean"

    sup_quals = [e.qual for e in supporting if e.qual is not None]
    if (
        strand_status == "both"
        and trace_quality == "clean"
        and sup_quals
        and min(sup_quals) >= cfg.high_qual
    ):
        confidence = "high"
    elif (
        strand_status in ("one-read", "single")
        and sup_quals
        and max(sup_quals) < cfg.low_qual
    ):
        confidence = "low"
    else:
        confidence = "review"

    return Variant(
        id=f"{sample.key}:{key.kind}:{key.ref_pos}:{key.ref_bases}>{key.alt_bases}",
        kind=key.kind,
        ref_pos=key.ref_pos,
        ref_bases=key.ref_bases,
        alt_bases=key.alt_bases,
        mixed=mixed,
        alt_fraction=alt_fraction,
        evidence=evidence,
        strand_status=strand_status,
        trace_quality=trace_quality,
        confidence=confidence,
    )
