import numpy as np

from sanger_seek.core.trim import mott_trim


def test_trims_low_quality_ends():
    quals = np.array([5] * 20 + [45] * 100 + [6] * 20, dtype=np.uint8)
    start, end = mott_trim(quals, len(quals))
    assert start == 20
    assert end == 120


def test_all_bad_returns_empty():
    quals = np.array([4] * 50, dtype=np.uint8)
    assert mott_trim(quals, 50) == (0, 0)


def test_no_quals_keeps_everything():
    assert mott_trim(None, 42) == (0, 42)


def test_interior_dip_is_kept():
    quals = np.array([5] * 10 + [45] * 50 + [12] * 5 + [45] * 50 + [5] * 10, dtype=np.uint8)
    start, end = mott_trim(quals, len(quals))
    assert start == 10
    assert end == 115
