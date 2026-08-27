"""SEQ parsing: plain base calls and Mutation Surveyor annotated references."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .dna import VALID_BASES


@dataclass
class MutationSurveyorSeq:
    """Structured content of a Mutation Surveyor annotated ``.seq`` file."""

    metadata: dict[str, str] = field(default_factory=dict)
    sequence: str = ""

    @property
    def gene(self) -> Optional[str]:
        return self.metadata.get("Gene")

    @property
    def translation(self) -> Optional[str]:
        return self.metadata.get("Translation")


METADATA_RE = re.compile(r"^\s*/(?P<key>.*?)\s*=\s*(?P<value>.*)$")
SEQUENCE_RE = re.compile(
    r"^\s*(?P<position>\d+)\s+(?P<bases>[ACGTUNRYKMSWBDHVX\-\s]+)\s*$",
    re.IGNORECASE,
)


def looks_like_mutation_surveyor(text: str) -> bool:
    """Require both slash metadata and a numbered nucleotide row."""
    lines = text.splitlines()
    return any(METADATA_RE.match(line) for line in lines) and any(
        SEQUENCE_RE.match(line) for line in lines
    )


def parse_mutation_surveyor_text(text: str, name: str = "") -> MutationSurveyorSeq:
    """Parse Mutation Surveyor metadata and numbered nucleotide rows.

    Multiline quoted values are consumed as one metadata field, preventing
    numbered protein-translation rows from being mistaken for DNA.
    """
    lines = text.splitlines()
    metadata: dict[str, str] = {}
    sequence_parts: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        metadata_match = METADATA_RE.match(line)
        if metadata_match:
            key = metadata_match.group("key").strip()
            value = metadata_match.group("value").strip()
            if value.startswith('"') and not _quoted_value_is_complete(value):
                collected = [value]
                i += 1
                while i < len(lines):
                    collected.append(lines[i])
                    if _quoted_value_is_complete("\n".join(collected)):
                        break
                    i += 1
                value = "\n".join(collected)
            metadata[key] = _clean_metadata_value(value)
            i += 1
            continue

        sequence_match = SEQUENCE_RE.match(line)
        if sequence_match:
            raw = re.sub(r"\s+", "", sequence_match.group("bases")).upper()
            # Reference coordinates must remain one base per source symbol.
            sequence_parts.append(raw.replace("U", "T").replace("X", "N").replace("-", "N"))
        i += 1

    sequence = "".join(sequence_parts)
    if not sequence:
        raise ValueError(f"{name or 'Mutation Surveyor SEQ'}: no numbered DNA sequence found")
    bad = sorted(set(sequence) - VALID_BASES)
    if bad:
        raise ValueError(
            f"{name or 'Mutation Surveyor SEQ'}: unexpected nucleotide symbols: "
            f"{''.join(bad)[:10]}"
        )
    return MutationSurveyorSeq(metadata=metadata, sequence=sequence)


def parse_mutation_surveyor_seq(path: str | Path) -> MutationSurveyorSeq:
    p = Path(path)
    return parse_mutation_surveyor_text(p.read_text(errors="replace"), name=p.name)


def _quoted_value_is_complete(value: str) -> bool:
    quote_count = 0
    escaped = False
    for char in value.rstrip():
        if char == "\\" and not escaped:
            escaped = True
            continue
        if char == '"' and not escaped:
            quote_count += 1
        escaped = False
    return quote_count >= 2 and quote_count % 2 == 0


def _clean_metadata_value(value: str) -> str:
    value = value.strip()
    if value.endswith(";"):
        value = value[:-1].rstrip()
    if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return value


def parse_seq_text(text: str, name: str = "") -> str:
    if looks_like_mutation_surveyor(text):
        return parse_mutation_surveyor_text(text, name=name).sequence
    lines = text.splitlines()
    body: list[str] = []
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith((">", ";", "#")):
            continue
        body.append(s)
    seq = "".join(body).upper()
    # Tolerate GCG/numbered formats: strip digits, whitespace and separators.
    seq = "".join(c for c in seq if not (c.isdigit() or c.isspace() or c in ".*"))
    bad = {c for c in seq if c not in VALID_BASES}
    if bad:
        if len(bad) > 4 or (len(seq) and sum(seq.count(c) for c in bad) > 0.1 * len(seq)):
            raise ValueError(
                f"{name or 'sequence'}: does not look like a nucleotide sequence "
                f"(unexpected characters: {''.join(sorted(bad))[:10]})"
            )
        seq = "".join(c if c in VALID_BASES else "N" for c in seq)
    if not seq:
        raise ValueError(f"{name or 'sequence'}: empty sequence")
    return seq


def load_seq(path: str | Path) -> str:
    p = Path(path)
    return parse_seq_text(p.read_text(errors="replace"), name=p.name)
