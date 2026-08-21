import numpy as np

from sanger_seek.core.abif import load_ab1
from sanger_seek.devtools.abif_writer import (
    tag_char,
    tag_pstring,
    tag_short,
    write_abif,
)


def test_roundtrip_minimal(tmp_path):
    calls = "ACGTN"
    ploc = [10, 22, 34, 46, 58]
    quals = bytes([40, 50, 12, 33, 2])
    data = {b: np.arange(70) + i * 100 for i, b in enumerate("GATC")}
    tags = [tag_char("FWO_", 1, "GATC"), tag_pstring("SMPL", 1, "testsample")]
    for i, b in enumerate("GATC"):
        tags.append(tag_short("DATA", 9 + i, data[b]))
    for num in (1, 2):
        tags.append(tag_char("PBAS", num, calls))
        tags.append(tag_short("PLOC", num, ploc))
        tags.append(tag_char("PCON", num, quals))

    path = tmp_path / "mini.ab1"
    write_abif(path, tags)

    trace = load_ab1(path)
    assert trace.calls == calls
    assert list(trace.ploc) == ploc
    assert list(trace.quals) == list(quals)
    assert trace.metadata.get("sample") == "testsample"
    for b in "GATC":
        assert np.array_equal(trace.channels[b], data[b])
    assert trace.n_samples == 70


def test_demo_ab1_loads(demo_dir):
    trace = load_ab1(demo_dir / "Sample001_F.ab1")
    n = len(trace.calls)
    assert n > 1000
    assert len(trace.ploc) == n and len(trace.quals) == n
    assert set(trace.channels) == set("ACGT")
    # peak locations must be within the trace
    assert trace.ploc.max() < trace.n_samples
