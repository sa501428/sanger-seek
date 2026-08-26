"""Read <-> reference alignment.

Primary engine is edlib (fast C++ infix alignment, IUPAC-aware via
additionalEqualities). A pure-Python/NumPy fallback implements the same
unit-cost infix ("HW") alignment so the app still works without edlib.

Conventions (SAM-like, verified against edlib output):
    '=' match        consumes read + ref
    'X' mismatch     consumes read + ref
    'I' insertion    consumes read only (extra bases in read)
    'D' deletion     consumes ref only (bases missing from read)
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from .dna import AMBIG_EQUALITIES, bases_match, expand, revcomp
from .model import Config, Read, ReadAlignment

try:
    import edlib  # type: ignore

    HAVE_EDLIB = True
except ImportError:  # pragma: no cover - depends on environment
    edlib = None
    HAVE_EDLIB = False

_CIGAR_RE = re.compile(r"(\d+)([=XIDM])")

MIN_ALIGN_LEN = 20


@dataclass
class RawAlignment:
    edit_distance: int
    ref_start: int
    ref_end: int                      # exclusive
    ops: list[tuple[str, int]]        # '=', 'X', 'I', 'D' runs (no 'M')


def parse_cigar(cigar: str) -> list[tuple[str, int]]:
    return [(op, int(n)) for n, op in _CIGAR_RE.findall(cigar)]


def _resolve_m_ops(
    ops: list[tuple[str, int]], query: str, target: str, ref_start: int
) -> list[tuple[str, int]]:
    """Split 'M' runs into '='/'X' and merge adjacent equal ops."""
    out: list[tuple[str, int]] = []

    def push(op: str, n: int) -> None:
        if n <= 0:
            return
        if out and out[-1][0] == op:
            out[-1] = (op, out[-1][1] + n)
        else:
            out.append((op, n))

    qi, ti = 0, ref_start
    for op, n in ops:
        if op == "M":
            for k in range(n):
                push("=" if bases_match(query[qi + k], target[ti + k]) else "X", 1)
            qi += n
            ti += n
        else:
            push(op, n)
            if op in ("=", "X"):
                qi += n
                ti += n
            elif op == "I":
                qi += n
            else:
                ti += n
    return out


# ---------------------------------------------------------------------------
# Pure-Python fallback (unit-cost infix alignment with IUPAC awareness)
# ---------------------------------------------------------------------------

def _match_masks(query: str, target: str) -> dict[str, np.ndarray]:
    tarr = np.frombuffer(target.encode("ascii", "replace"), dtype="S1")
    masks: dict[str, np.ndarray] = {}
    for ch in set(query):
        eq = [b.encode() for b in expand(ch)] or [ch.encode()]
        mask = np.zeros(len(target), dtype=bool)
        for t in set(target):
            if bases_match(ch, t):
                mask |= tarr == t.encode()
        masks[ch] = mask
        del eq
    return masks


def _window_target(query: str, target: str, k: int = 12) -> tuple[int, str]:
    """K-mer seed to restrict very long targets to a window around the read."""
    if len(target) <= 20000:
        return 0, target
    index: dict[str, int] = {}
    for j in range(len(target) - k + 1):
        index.setdefault(target[j : j + k], j)
    diags: list[int] = []
    for i in range(0, len(query) - k + 1, 4):
        j = index.get(query[i : i + k])
        if j is not None:
            diags.append(j - i)
    if not diags:
        return 0, target
    diags.sort()
    d = diags[len(diags) // 2]
    lo = max(0, d - 500)
    hi = min(len(target), d + len(query) + 500)
    return lo, target[lo:hi]


def _py_align(query: str, target: str) -> RawAlignment | None:
    offset, target = _window_target(query, target)
    m, n = len(query), len(target)
    if m == 0 or n == 0:
        return None
    masks = _match_masks(query, target)
    D = np.zeros((m + 1, n + 1), dtype=np.int32)
    D[0, :] = 0
    D[:, 0] = np.arange(m + 1)
    steps = np.arange(n + 1, dtype=np.int32)
    for i in range(1, m + 1):
        prev = D[i - 1]
        cost = np.where(masks[query[i - 1]], 0, 1).astype(np.int32)
        tmp = np.minimum(prev[:-1] + cost, prev[1:] + 1)
        arr = np.concatenate(([np.int32(i)], tmp))
        D[i] = np.minimum.accumulate(arr - steps) + steps

    j = int(np.argmin(D[m]))
    dist = int(D[m, j])
    ops_rev: list[str] = []
    i = m
    while i > 0:
        if j > 0:
            diag_cost = 0 if bases_match(query[i - 1], target[j - 1]) else 1
            if D[i, j] == D[i - 1, j - 1] + diag_cost:
                ops_rev.append("=" if diag_cost == 0 else "X")
                i -= 1
                j -= 1
                continue
        if D[i, j] == D[i - 1, j] + 1:
            ops_rev.append("I")
            i -= 1
            continue
        ops_rev.append("D")
        j -= 1

    runs: list[tuple[str, int]] = []
    for op in reversed(ops_rev):
        if runs and runs[-1][0] == op:
            runs[-1] = (op, runs[-1][1] + 1)
        else:
            runs.append((op, 1))
    ref_start = offset + j
    span = sum(n_ for op, n_ in runs if op in ("=", "X", "D"))
    return RawAlignment(dist, ref_start, ref_start + span, runs)


def _py_distance(query: str, target: str) -> int:
    offset, target = _window_target(query, target)
    del offset
    m, n = len(query), len(target)
    masks = _match_masks(query, target)
    prev = np.zeros(n + 1, dtype=np.int32)
    steps = np.arange(n + 1, dtype=np.int32)
    for i in range(1, m + 1):
        cost = np.where(masks[query[i - 1]], 0, 1).astype(np.int32)
        tmp = np.minimum(prev[:-1] + cost, prev[1:] + 1)
        arr = np.concatenate(([np.int32(i)], tmp))
        prev = np.minimum.accumulate(arr - steps) + steps
    return int(prev.min())


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def infix_distance(query: str, target: str) -> int:
    """Edit distance of the best infix alignment of query inside target."""
    if HAVE_EDLIB:
        r = edlib.align(
            query, target, mode="HW", task="distance",
            additionalEqualities=AMBIG_EQUALITIES,
        )
        return int(r["editDistance"])
    return _py_distance(query, target)


def global_distance(a: str, b: str) -> int:
    """Global (NW) edit distance between two sequences."""
    if HAVE_EDLIB:
        return int(edlib.align(a, b, mode="NW", task="distance")["editDistance"])
    if len(a) == len(b):
        return sum(1 for x, y in zip(a, b) if x != y)
    return _py_distance(a, b) if len(a) <= len(b) else _py_distance(b, a)


def align_infix(query: str, target: str) -> RawAlignment | None:
    """Best infix alignment of query inside target, with full op path."""
    if not query or not target:
        return None
    if HAVE_EDLIB:
        r = edlib.align(
            query, target, mode="HW", task="path",
            additionalEqualities=AMBIG_EQUALITIES,
        )
        if r["editDistance"] < 0 or not r.get("cigar"):
            return None
        tstart, tend = r["locations"][0]
        ops = parse_cigar(r["cigar"])
        ops = _resolve_m_ops(ops, query, target, int(tstart))
        return RawAlignment(int(r["editDistance"]), int(tstart), int(tend) + 1, ops)
    return _py_align(query, target)


def orient_and_align(read: Read, ref_seq: str, cfg: Config) -> ReadAlignment | None:
    """Decide read orientation and align its trimmed window to the reference.

    Returns a ReadAlignment in reference orientation, or None when the read
    does not align acceptably in either orientation.
    """
    n = read.n
    ts, te = read.trim
    if te - ts < MIN_ALIGN_LEN:
        return None

    fwd = read.calls
    rev = revcomp(read.calls)
    # Trimmed window expressed in each orientation's coordinates.
    windows = {
        "F": (fwd[ts:te], ts),
        "R": (rev[n - te : n - ts], n - te),
    }

    dists = {o: infix_distance(q, ref_seq) for o, (q, _) in windows.items()}
    order = (
        [read.orientation_override]
        if read.orientation_override in dists
        else sorted(dists, key=lambda o: dists[o])
    )
    if len(order) > 1 and dists[order[0]] == dists[order[1]] and read.orientation_hint in dists:
        order = [read.orientation_hint] + [o for o in order if o != read.orientation_hint]

    for orientation in order:
        query, win_start = windows[orientation]
        raw, query_offset, aligned_query_len = _align_read_window(query, ref_seq)
        if raw is None:
            continue
        if raw.edit_distance / max(len(query), 1) > cfg.max_align_frac_dist:
            continue
        return _build_alignment(
            raw,
            orientation,
            win_start + query_offset,
            aligned_query_len,
            n,
            len(ref_seq),
        )
    return None


def _align_read_window(
    query: str, ref_seq: str
) -> tuple[RawAlignment | None, int, int]:
    """Align either sequence as the contained interval.

    Amplicon references are often shorter than the chromatogram because the
    read includes primer or flanking sequence. In that case, align the
    reference inside the read, invert the operation path, and clip the read
    overhangs instead of reporting them as giant insertions.
    """
    if len(query) <= len(ref_seq):
        return align_infix(query, ref_seq), 0, len(query)
    swapped = align_infix(ref_seq, query)
    if swapped is None:
        return None, 0, 0
    inverse = [
        ("D" if op == "I" else "I" if op == "D" else op, length)
        for op, length in swapped.ops
    ]
    raw = RawAlignment(
        edit_distance=swapped.edit_distance,
        ref_start=0,
        ref_end=len(ref_seq),
        ops=inverse,
    )
    query_offset = swapped.ref_start
    aligned_query_len = swapped.ref_end - swapped.ref_start
    return raw, query_offset, aligned_query_len


def _build_alignment(
    raw: RawAlignment,
    orientation: str,
    win_start: int,
    win_len: int,
    read_len: int,
    ref_len: int,
) -> ReadAlignment:
    read_to_ref = np.full(read_len, -1, dtype=np.int64)
    ref_to_read = np.full(raw.ref_end - raw.ref_start, -1, dtype=np.int64)

    rc = win_start          # oriented read cursor
    fc = raw.ref_start      # reference cursor
    matches = 0
    columns = 0
    for op, k in raw.ops:
        if op in ("=", "X"):
            for t in range(k):
                read_to_ref[rc + t] = fc + t
                ref_to_read[fc + t - raw.ref_start] = rc + t
            if op == "=":
                matches += k
            rc += k
            fc += k
            columns += k
        elif op == "I":
            rc += k
            columns += k
        elif op == "D":
            fc += k
            columns += k

    return ReadAlignment(
        orientation=orientation,
        edit_distance=raw.edit_distance,
        identity=matches / columns if columns else 0.0,
        ref_start=raw.ref_start,
        ref_end=raw.ref_end,
        read_start=win_start,
        read_end=win_start + win_len,
        ops=raw.ops,
        read_to_ref=read_to_ref,
        ref_to_read=ref_to_read,
    )
