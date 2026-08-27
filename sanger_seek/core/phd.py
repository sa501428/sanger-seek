"""Parser for optional Phred ``.phd`` / ``.phd.1`` base-call files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class PhdData:
    calls: str
    quals: np.ndarray
    ploc: np.ndarray


def parse_phd_text(text: str) -> PhdData:
    """Read the base, quality and peak-position columns in ``BEGIN_DNA``.

    Phred permits extra columns after the three standard fields, so they are
    deliberately ignored.  Calls are normalized to uppercase.
    """
    in_dna = False
    calls: list[str] = []
    quals: list[int] = []
    ploc: list[int] = []
    for raw in text.splitlines():
        line = raw.strip()
        if line == "BEGIN_DNA":
            in_dna = True
            continue
        if line == "END_DNA":
            break
        if not in_dna or not line:
            continue
        fields = line.split()
        if len(fields) < 3 or len(fields[0]) != 1 or fields[0].upper() not in "ACGTNRYKMSWBDHV":
            raise ValueError(f"invalid PHD DNA row: {line[:80]}")
        try:
            quality = int(fields[1])
            peak = int(fields[2])
        except ValueError as exc:
            raise ValueError(f"invalid PHD quality/peak row: {line[:80]}") from exc
        calls.append(fields[0].upper())
        quals.append(max(0, min(quality, 255)))
        ploc.append(max(0, peak))
    if not calls:
        raise ValueError("no base calls found between BEGIN_DNA and END_DNA")
    return PhdData(
        calls="".join(calls),
        quals=np.asarray(quals, dtype=np.uint8),
        ploc=np.asarray(ploc, dtype=np.int32),
    )


def load_phd(path: str | Path) -> PhdData:
    p = Path(path)
    try:
        return parse_phd_text(p.read_text(errors="replace"))
    except Exception as exc:
        raise ValueError(f"{p.name}: {exc}") from exc
