"""Reference loading: FASTA (sequence only) or GenBank (sequence + CDS features)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio import SeqIO

from .dna import complement
from .model import CDSFeature, Reference

GENBANK_EXTS = {".gb", ".gbk", ".genbank", ".gbff"}
FASTA_EXTS = {".fa", ".fasta", ".fna", ".ffn", ".txt"}


def _qual1(feature, key: str, default: str = "") -> str:
    v = feature.qualifiers.get(key)
    return str(v[0]) if v else default


def _cds_from_feature(feature, genome: str) -> CDSFeature:
    parts: list[tuple[int, int]] = []
    order: list[np.ndarray] = []
    for part in feature.location.parts:
        s, e = int(part.start), int(part.end)
        parts.append((s, e))
        idx = np.arange(s, e, dtype=np.int64)
        if part.strand == -1:
            idx = idx[::-1]
        order.append(idx)
    genomic_order = np.concatenate(order) if order else np.zeros(0, dtype=np.int64)

    strand = -1 if feature.location.strand == -1 else 1
    cds_seq = "".join(
        complement(genome[g]) if strand == -1 else genome[g] for g in genomic_order
    )

    codon_start = 1
    try:
        codon_start = int(_qual1(feature, "codon_start", "1"))
    except ValueError:
        pass

    gene = (
        _qual1(feature, "gene")
        or _qual1(feature, "locus_tag")
        or _qual1(feature, "product")
        or "CDS"
    )
    return CDSFeature(
        gene=gene,
        product=_qual1(feature, "product"),
        strand=strand,
        codon_start=min(max(codon_start, 1), 3),
        parts=parts,
        genomic_order=genomic_order,
        cds_seq=cds_seq,
    )


def load_reference(path: str | Path) -> Reference:
    p = Path(path)
    fmt = "genbank" if p.suffix.lower() in GENBANK_EXTS else "fasta"
    records = list(SeqIO.parse(str(p), fmt))
    if not records:
        raise ValueError(f"{p.name}: no sequence records found")
    record = records[0]
    seq = str(record.seq).upper()
    if not seq:
        raise ValueError(f"{p.name}: empty reference sequence")

    cds: list[CDSFeature] = []
    if fmt == "genbank":
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    cds.append(_cds_from_feature(feature, seq))
                except Exception:
                    continue

    name = record.id if record.id and record.id != "<unknown id>" else p.stem
    return Reference(
        name=name,
        seq=seq,
        path=str(p),
        source=fmt,
        cds=cds,
        description=getattr(record, "description", "") or "",
    )
