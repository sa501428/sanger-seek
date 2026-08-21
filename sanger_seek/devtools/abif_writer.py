"""Minimal ABIF (.ab1) writer — for generating synthetic demo/test traces.

Implements the ABIF container format: 128-byte header with a root directory
entry, data blocks, then the tag directory. Data <= 4 bytes is stored inline
in the entry's dataoffset field, as per the spec.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

# element type codes
ET_BYTE = 1
ET_CHAR = 2
ET_SHORT = 4
ET_LONG = 5
ET_PSTRING = 18
ET_CSTRING = 19


@dataclass
class Tag:
    name: str          # 4 chars
    number: int
    etype: int
    esize: int
    nelem: int
    data: bytes


def tag_char(name: str, number: int, value: str | bytes) -> Tag:
    b = value.encode("ascii") if isinstance(value, str) else bytes(value)
    return Tag(name, number, ET_CHAR, 1, len(b), b)


def tag_byte(name: str, number: int, values) -> Tag:
    b = bytes(bytearray(int(v) & 0xFF for v in values))
    return Tag(name, number, ET_BYTE, 1, len(b), b)


def tag_short(name: str, number: int, values) -> Tag:
    vals = [int(v) for v in values]
    return Tag(name, number, ET_SHORT, 2, len(vals), struct.pack(f">{len(vals)}h", *vals))


def tag_long(name: str, number: int, values) -> Tag:
    vals = [int(v) for v in values]
    return Tag(name, number, ET_LONG, 4, len(vals), struct.pack(f">{len(vals)}i", *vals))


def tag_pstring(name: str, number: int, value: str) -> Tag:
    b = value.encode("ascii")[:255]
    return Tag(name, number, ET_PSTRING, 1, len(b) + 1, bytes([len(b)]) + b)


def write_abif(path: str | Path, tags: list[Tag]) -> None:
    tags = sorted(tags, key=lambda t: (t.name, t.number))
    blocks: list[bytes] = []
    offsets: list[int] = []
    pos = 128  # data starts after the header
    for t in tags:
        if len(t.data) > 4:
            offsets.append(pos)
            blocks.append(t.data)
            pos += len(t.data)
        else:
            offsets.append(-1)  # inline
    dir_offset = pos

    entries = bytearray()
    for t, off in zip(tags, offsets):
        name = t.name.encode("ascii").ljust(4)[:4]
        if off < 0:
            inline = t.data.ljust(4, b"\x00")
            entries += struct.pack(">4sihhii4si", name, t.number, t.etype, t.esize,
                                   t.nelem, len(t.data), inline, 0)
        else:
            entries += struct.pack(">4sihhiiii", name, t.number, t.etype, t.esize,
                                   t.nelem, len(t.data), off, 0)

    header = bytearray(128)
    header[0:4] = b"ABIF"
    struct.pack_into(">H", header, 4, 101)  # version 1.01
    struct.pack_into(
        ">4sihhiiii", header, 6,
        b"tdir", 1, 1023, 28, len(tags), 28 * len(tags), dir_offset, 0,
    )

    with open(path, "wb") as fh:
        fh.write(header)
        for b in blocks:
            fh.write(b)
        fh.write(entries)
