""".seq file parsing: plain-text exported base calls, optionally FASTA-headered."""

from __future__ import annotations

from pathlib import Path

from .dna import VALID_BASES


def parse_seq_text(text: str, name: str = "") -> str:
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
