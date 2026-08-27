"""CSV export of the variant review table."""

from __future__ import annotations

import csv
from pathlib import Path

from .model import Reference, Sample, Variant

HEADER = [
    "case", "variant", "protein", "effect", "gene", "type",
    "ref_pos_1based", "ref", "alt", "mixed", "alt_fraction",
    "fwd", "rev", "strand_status", "trace_quality", "confidence",
    "evidence",
]


def _strand_mark(v: Variant, orientation: str) -> str:
    evs = [e for e in v.evidence if e.orientation == orientation]
    if not any(e.covered for e in evs):
        return "-"
    if any(e.supports for e in evs):
        return "yes"
    if all(e.supports is None for e in evs if e.covered):
        return "?"
    return "no"


def _evidence_text(v: Variant) -> str:
    parts = []
    for e in v.evidence:
        if not e.covered:
            parts.append(f"{e.read_label}: not covered")
            continue
        bits = [f"{e.read_label}: {e.call or '?'}"]
        if e.qual is not None:
            bits.append(f"Q{e.qual}")
        if e.ratio is not None:
            bits.append(f"2nd/1st {e.ratio:.2f}")
        bits.append("supports" if e.supports else "no-support" if e.supports is False else "ambiguous")
        parts.append(" ".join(bits))
    return "; ".join(parts)


def variant_row(sample: Sample, v: Variant) -> list[str]:
    return [
        sample.name,
        v.label,
        v.protein or "",
        v.effect or "",
        v.gene or "",
        v.kind,
        str(v.ref_pos + 1),
        v.ref_bases,
        v.alt_bases,
        "yes" if v.mixed else "no",
        f"{v.alt_fraction:.2f}" if v.alt_fraction is not None else "",
        _strand_mark(v, "F"),
        _strand_mark(v, "R"),
        v.strand_status,
        v.trace_quality,
        v.confidence,
        _evidence_text(v),
    ]


def export_variants_csv(
    path: str | Path,
    samples: list[Sample],
    reference: Reference | None = None,
    variants_filter=None,
) -> int:
    """Write variants for the given samples; returns row count."""
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        if reference is not None:
            w.writerow([f"# reference: {reference.name} ({reference.path})"])
        w.writerow(HEADER)
        for s in samples:
            for v in s.variants:
                if variants_filter is not None and not variants_filter(v):
                    continue
                w.writerow(variant_row(s, v))
                n += 1
    return n
