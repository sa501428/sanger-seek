"""Quality trimming using the Mott algorithm (as in phred/trace tuner).

Each base contributes (cutoff - error_prob); the trim window is the
maximum-sum contiguous segment (Kadane). Returns [start, end) on the
original call order; (0, 0) when nothing passes.
"""

from __future__ import annotations

import numpy as np


def mott_trim(quals: np.ndarray | None, n: int, cutoff: float = 0.05) -> tuple[int, int]:
    if quals is None or len(quals) == 0:
        return (0, n)
    q = np.asarray(quals, dtype=np.float64)[:n]
    scores = cutoff - np.power(10.0, q / -10.0)

    best_sum = 0.0
    best = (0, 0)
    cur_sum = 0.0
    cur_start = 0
    for i, s in enumerate(scores):
        cur_sum += s
        if cur_sum <= 0.0:
            cur_sum = 0.0
            cur_start = i + 1
        elif cur_sum > best_sum:
            best_sum = cur_sum
            best = (cur_start, i + 1)
    return best
