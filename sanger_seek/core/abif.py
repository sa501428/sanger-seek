"""AB1 (ABIF) loading via Biopython.

The AB1 file is the authoritative source for the electropherogram. We extract
the analyzed traces (DATA9-12 in FWO_ channel order), base calls (PBAS),
peak locations (PLOC), per-base qualities (PCON) and light metadata.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio import SeqIO

from .model import TraceData

_META_TAGS = {
    "SMPL1": "sample",
    "MCHN1": "instrument",
    "MODL1": "model",
    "DySN1": "dye_set",
    "RMdN1": "run_module",
    "TUBE1": "well",
    "CMNT1": "comment",
}


def _text(v) -> str:
    if isinstance(v, (bytes, bytearray)):
        return v.decode("ascii", errors="replace")
    return str(v)


def _first(raw: dict, *keys):
    for k in keys:
        if k in raw and raw[k] is not None:
            return raw[k]
    return None


def load_ab1(path: str | Path) -> TraceData:
    record = SeqIO.read(str(path), "abi")
    raw = record.annotations.get("abif_raw")
    if not raw:
        raise ValueError(f"{path}: no ABIF tag data found")

    order = _text(raw.get("FWO_1", b"GATC")).upper()
    channels: dict[str, np.ndarray] = {}
    for i, base in enumerate(order[:4]):
        data = raw.get(f"DATA{9 + i}")
        if data is None:
            raise ValueError(f"{path}: missing analyzed trace DATA{9 + i}")
        channels[base] = np.asarray(data, dtype=np.int32)

    calls_raw = _first(raw, "PBAS2", "PBAS1")
    ploc_raw = _first(raw, "PLOC2", "PLOC1")
    pcon_raw = _first(raw, "PCON2", "PCON1")
    if calls_raw is None or ploc_raw is None:
        raise ValueError(f"{path}: missing base calls (PBAS) or peak locations (PLOC)")

    calls = _text(calls_raw).upper()
    ploc = np.asarray(ploc_raw, dtype=np.int64)
    if pcon_raw is None:
        quals = np.zeros(len(calls), dtype=np.uint8)
    elif isinstance(pcon_raw, (bytes, bytearray)):
        quals = np.frombuffer(bytes(pcon_raw), dtype=np.uint8).copy()
    else:
        quals = np.asarray(pcon_raw, dtype=np.uint8)

    n = min(len(calls), len(ploc), len(quals))
    calls, ploc, quals = calls[:n], ploc[:n], quals[:n]

    metadata: dict[str, str] = {}
    for tag, name in _META_TAGS.items():
        if tag in raw:
            try:
                metadata[name] = _text(raw[tag]).strip()
            except Exception:
                pass

    return TraceData(
        channels=channels, calls=calls, quals=quals, ploc=ploc, metadata=metadata
    )
