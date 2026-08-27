"""File classification and sample/read pairing.

Files are grouped into samples by basename after stripping a trailing
direction token (_F, -R, _FWD, _rev, ...). `sample_F.seq` pairs with
`sample_F.ab1` (same stem = same read). Direction tokens are only a hint;
final orientation is decided by alignment.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

AB1_EXTS = {".ab1", ".abi", ".ab"}
SEQ_EXTS = {".seq"}
PHD_RE = re.compile(r"\.phd(?:\.\d+)?$", re.IGNORECASE)
FASTA_EXTS = {".fa", ".fasta", ".fna", ".ffn"}
GENBANK_EXTS = {".gb", ".gbk", ".genbank", ".gbff"}

_FWD_TOKENS = {"f", "fw", "fwd", "for", "forward"}
_REV_TOKENS = {"r", "rv", "rev", "reverse"}
_DIR_RE = re.compile(
    r"[._\- ](f|r|fw|rv|fwd|rev|for|forward|reverse)(\d{0,2})$", re.IGNORECASE
)


def classify(path: str | Path) -> str | None:
    p = Path(path)
    if PHD_RE.search(p.name):
        return "phd"
    ext = p.suffix.lower()
    if ext in AB1_EXTS:
        return "ab1"
    if ext in SEQ_EXTS:
        return "seq"
    if ext in FASTA_EXTS:
        return "fasta"
    if ext in GENBANK_EXTS:
        return "genbank"
    return None


def split_direction(stem: str) -> tuple[str, str | None]:
    """Return (sample_key, orientation_hint) for a file stem."""
    m = _DIR_RE.search(stem)
    if not m:
        return stem, None
    token = m.group(1).lower()
    hint = "F" if token in _FWD_TOKENS else "R" if token in _REV_TOKENS else None
    return stem[: m.start()], hint


@dataclass
class ReadFiles:
    stem: str
    ab1: Path | None = None
    seq: Path | None = None
    phd: Path | None = None
    hint: str | None = None


@dataclass
class ScanResult:
    samples: dict[str, dict[str, ReadFiles]] = field(default_factory=dict)
    references: list[Path] = field(default_factory=list)
    ignored: list[Path] = field(default_factory=list)

    def merge(self, other: "ScanResult") -> None:
        for skey, reads in other.samples.items():
            mine = self.samples.setdefault(skey, {})
            for rkey, rf in reads.items():
                if rkey in mine:
                    mine[rkey].ab1 = mine[rkey].ab1 or rf.ab1
                    mine[rkey].seq = mine[rkey].seq or rf.seq
                    mine[rkey].phd = mine[rkey].phd or rf.phd
                    mine[rkey].hint = mine[rkey].hint or rf.hint
                else:
                    mine[rkey] = rf
        self.references.extend(other.references)
        self.ignored.extend(other.ignored)


def scan_paths(paths: list[str | Path]) -> ScanResult:
    result = ScanResult()
    files: list[Path] = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            files.extend(sorted(q for q in p.iterdir() if q.is_file()))
        else:
            files.append(p)

    for f in files:
        kind = classify(f)
        if kind is None:
            result.ignored.append(f)
            continue
        if kind in ("fasta", "genbank"):
            result.references.append(f)
            continue
        if kind == "phd":
            stem = PHD_RE.sub("", f.name)
            # Common Phred output is ``read.ab1.phd.1``.
            if Path(stem).suffix.lower() in AB1_EXTS:
                stem = Path(stem).stem
        else:
            stem = f.stem
        sample_key, hint = split_direction(stem)
        rkey = stem.lower()
        reads = result.samples.setdefault(sample_key, {})
        rf = reads.setdefault(rkey, ReadFiles(stem=stem, hint=hint))
        if kind == "ab1":
            rf.ab1 = f
        elif kind == "seq":
            rf.seq = f
        else:
            rf.phd = f
        if rf.hint is None:
            rf.hint = hint
    return result
