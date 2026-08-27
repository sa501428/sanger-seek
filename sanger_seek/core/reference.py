"""Reference loading: plain/Mutation Surveyor SEQ, FASTA, or GenBank."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from Bio import SeqIO

from .dna import complement
from .model import CDSFeature, Reference
from .seqfile import (
    looks_like_mutation_surveyor,
    load_seq,
    parse_mutation_surveyor_text,
)

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
    if p.suffix.lower() == ".seq":
        # Some GenBank exports retain the generic .seq extension. Sniff the
        # record header before treating it as plain base-call text.
        text = p.read_text(errors="replace")
        if not text.lstrip().startswith("LOCUS"):
            if looks_like_mutation_surveyor(text):
                parsed = parse_mutation_surveyor_text(text, name=p.name)
                cds = _mutation_surveyor_cds(parsed.metadata, parsed.sequence, p.stem)
                description = (
                    parsed.metadata.get("Exon_And_Note")
                    or parsed.metadata.get("Exon And Note")
                    or "Mutation Surveyor annotated reference"
                )
                return Reference(
                    name=parsed.gene or p.stem,
                    seq=parsed.sequence,
                    path=str(p),
                    source="mutation-surveyor",
                    cds=[cds] if cds else [],
                    description=description,
                    metadata=parsed.metadata,
                )
            seq = load_seq(p)
            if not seq:
                raise ValueError(f"{p.name}: empty reference sequence")
            return Reference(name=p.stem, seq=seq, path=str(p), source="seq")
        fmt = "genbank"
    else:
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


def _mutation_surveyor_cds(
    metadata: dict[str, str], genome: str, fallback_gene: str
) -> CDSFeature | None:
    """Convert simple Mutation Surveyor CDS coordinates to a CDS feature."""
    raw_cds = metadata.get("CDS", "")
    match = re.search(r"<?(\d+)\s*\.\.\s*>?(\d+)", raw_cds)
    if not match:
        return None
    start_1, end_1 = (int(value) for value in match.groups())
    start = max(start_1 - 1, 0)
    end = min(end_1, len(genome))
    if start >= end:
        return None
    strand = -1 if "complement" in raw_cds.lower() else 1
    if strand == 1:
        genomic_order = np.arange(start, end, dtype=np.int64)
        cds_seq = genome[start:end]
    else:
        genomic_order = np.arange(end - 1, start - 1, -1, dtype=np.int64)
        cds_seq = "".join(complement(genome[g]) for g in genomic_order)
    try:
        codon_start = int(metadata.get("Reading Frame (1,2,3)", "1"))
    except ValueError:
        codon_start = 1
    return CDSFeature(
        gene=metadata.get("Gene") or fallback_gene,
        product=metadata.get("Exon_And_Note", ""),
        strand=strand,
        codon_start=min(max(codon_start, 1), 3),
        parts=[(start, end)],
        genomic_order=genomic_order,
        cds_seq=cds_seq,
    )
