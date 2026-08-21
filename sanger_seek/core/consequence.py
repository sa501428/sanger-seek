"""Coding consequence annotation (best-effort HGVS-like c. / p. labels).

Requires a GenBank reference with CDS features. Handles multi-exon joins,
reverse-strand CDS and codon_start offsets. Labels are for review purposes,
not validated clinical nomenclature.
"""

from __future__ import annotations

from .dna import aa3, complement, revcomp, translate_codon
from .model import CDSFeature, Reference, Variant


def annotate_variants(variants: list[Variant], reference: Reference) -> None:
    if not reference.cds:
        return
    for v in variants:
        _annotate(v, reference)


def _feature_for(v: Variant, reference: Reference) -> CDSFeature | None:
    positions: list[int] = []
    if v.kind == "snv":
        positions = [v.ref_pos]
    elif v.kind == "del":
        positions = list(range(v.ref_pos, v.ref_pos + len(v.ref_bases)))
    else:
        positions = [v.ref_pos - 1, v.ref_pos]
    for p in positions:
        f = reference.feature_at(p)
        if f is not None:
            return f
    return None


def _annotate(v: Variant, reference: Reference) -> None:
    feature = _feature_for(v, reference)
    if feature is None:
        v.effect = "non-coding"
        return
    v.gene = feature.gene
    if v.kind == "snv":
        _annotate_snv(v, feature)
    elif v.kind == "del":
        _annotate_del(v, feature)
    else:
        _annotate_ins(v, feature)


def _coding(feature: CDSFeature) -> str:
    return feature.cds_seq[feature.codon_start - 1 :]


def _cpos(feature: CDSFeature, gpos: int) -> int | None:
    """Position within the coding sequence (0-based), or None."""
    ci = feature.cds_index(gpos)
    if ci is None:
        return None
    cpos = ci - (feature.codon_start - 1)
    return cpos if cpos >= 0 else None


def _annotate_snv(v: Variant, feature: CDSFeature) -> None:
    cpos = _cpos(feature, v.ref_pos)
    if cpos is None:
        v.effect = "non-coding"
        return
    coding = _coding(feature)
    ref_c = coding[cpos]
    alt_c = v.alt_bases if feature.strand == 1 else complement(v.alt_bases)
    v.cdna = f"c.{cpos + 1}{ref_c}>{alt_c}"

    codon_i = cpos // 3
    codon = coding[codon_i * 3 : codon_i * 3 + 3]
    if len(codon) < 3:
        v.effect = "coding (partial codon)"
        return
    within = cpos % 3
    alt_codon = codon[:within] + alt_c + codon[within + 1 :]
    aa_ref = translate_codon(codon)
    aa_alt = translate_codon(alt_codon)
    n = codon_i + 1
    if aa_alt == aa_ref:
        v.effect = "synonymous"
        v.protein = f"p.{aa3(aa_ref)}{n}="
    elif aa_alt == "*":
        v.effect = "nonsense"
        v.protein = f"p.{aa3(aa_ref)}{n}Ter"
    elif aa_ref == "*":
        v.effect = "stop-loss"
        v.protein = f"p.Ter{n}{aa3(aa_alt)}ext"
    else:
        v.effect = "missense"
        v.protein = f"p.{aa3(aa_ref)}{n}{aa3(aa_alt)}"


def _annotate_del(v: Variant, feature: CDSFeature) -> None:
    cposs = [
        c
        for p in range(v.ref_pos, v.ref_pos + len(v.ref_bases))
        if (c := _cpos(feature, p)) is not None
    ]
    if not cposs:
        v.effect = "non-coding"
        return
    coding = _coding(feature)
    c1, c2 = min(cposs) + 1, max(cposs) + 1
    if len(v.ref_bases) == 1:
        v.cdna = f"c.{c1}del{coding[c1 - 1]}"
    else:
        v.cdna = f"c.{c1}_{c2}del"
    n_coding = len(cposs)
    codon1 = (c1 - 1) // 3 + 1
    aa_first = translate_codon(coding[(codon1 - 1) * 3 : (codon1 - 1) * 3 + 3])
    if n_coding % 3 != 0:
        v.effect = "frameshift"
        v.protein = f"p.{aa3(aa_first)}{codon1}fs"
    else:
        v.effect = "in-frame deletion"
        codon2 = (c2 - 1) // 3 + 1
        if codon1 == codon2:
            v.protein = f"p.{aa3(aa_first)}{codon1}del"
        else:
            aa_last = translate_codon(coding[(codon2 - 1) * 3 : (codon2 - 1) * 3 + 3])
            v.protein = f"p.{aa3(aa_first)}{codon1}_{aa3(aa_last)}{codon2}del"


def _annotate_ins(v: Variant, feature: CDSFeature) -> None:
    left = _cpos(feature, v.ref_pos - 1)
    right = _cpos(feature, v.ref_pos)
    if left is None and right is None:
        v.effect = "non-coding"
        return
    bases = v.alt_bases if feature.strand == 1 else revcomp(v.alt_bases)
    anchors = sorted(c for c in (left, right) if c is not None)
    c1 = anchors[0] + 1
    v.cdna = f"c.{c1}_{c1 + 1}ins{bases}"
    coding = _coding(feature)
    codon1 = (c1 - 1) // 3 + 1
    aa_first = translate_codon(coding[(codon1 - 1) * 3 : (codon1 - 1) * 3 + 3])
    if len(v.alt_bases) % 3 != 0:
        v.effect = "frameshift"
        v.protein = f"p.{aa3(aa_first)}{codon1}fs"
    else:
        v.effect = "in-frame insertion"
        v.protein = f"p.{aa3(aa_first)}{codon1}ins{len(bases) // 3}aa"
