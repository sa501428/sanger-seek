"""Per-called-base signal metrics from the raw chromatogram.

For each called base we look in a window around its peak location (bounded
by the midpoints to neighboring peaks) and record the primary peak height,
the strongest other channel (secondary), their ratio, and the local
baseline noise. These drive mixed-peak detection and trace-quality flags.
"""

from __future__ import annotations

import numpy as np

from .model import PeakMetrics, TraceData

BASES = "ACGT"


def peak_metrics(trace: TraceData) -> list[PeakMetrics]:
    ploc = np.asarray(trace.ploc, dtype=np.int64)
    n = len(ploc)
    if n == 0:
        return []

    n_samples = trace.n_samples
    channels = {b: trace.channels.get(b, np.zeros(n_samples, dtype=np.int32)) for b in BASES}

    # Window boundaries: midpoints between adjacent peak locations.
    mids = (ploc[1:] + ploc[:-1]) // 2
    lo = np.concatenate(([max(ploc[0] - 6, 0)], mids))
    hi = np.concatenate((mids, [min(ploc[-1] + 7, n_samples)]))
    lo = np.minimum(lo, ploc)                      # window always includes the peak
    hi = np.maximum(hi, ploc + 1)

    out: list[PeakMetrics] = []
    for i in range(n):
        a, b = int(lo[i]), int(hi[i])
        heights = {}
        floors = []
        for base in BASES:
            win = channels[base][a:b]
            if len(win) == 0:
                heights[base] = 0.0
                floors.append(0.0)
            else:
                heights[base] = float(win.max())
                floors.append(float(win.min()))
        called = trace.calls[i].upper()
        if called in heights:
            primary = heights[called]
            others = {k: v for k, v in heights.items() if k != called}
        else:
            # Ambiguity/N call: primary = tallest channel.
            primary_base = max(heights, key=heights.get)
            primary = heights[primary_base]
            others = {k: v for k, v in heights.items() if k != primary_base}
        secondary_base = max(others, key=others.get)
        secondary = others[secondary_base]
        noise = float(np.median(floors))
        # Ignore "secondary" signal that is indistinguishable from baseline.
        eff_secondary = max(secondary - noise, 0.0)
        eff_primary = max(primary - noise, 1e-9)
        ratio = min(eff_secondary / eff_primary, 9.99) if primary > 0 else 0.0
        out.append(
            PeakMetrics(
                called=called,
                primary_h=primary,
                secondary_h=secondary,
                ratio=round(ratio, 4),
                secondary_base=secondary_base,
                noise=noise,
            )
        )
    return out
