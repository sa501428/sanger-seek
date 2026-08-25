"""Offscreen smoke coverage for the assembled desktop application."""

import os
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QSettings  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from sanger_seek.app import create_application  # noqa: E402
from sanger_seek.ui.main_window import MainWindow  # noqa: E402


def _wait_until(app, predicate, timeout=10.0):
    deadline = time.monotonic() + timeout
    while not predicate() and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(0.01)
    assert predicate(), "desktop operation timed out"


def test_desktop_demo_analysis_and_screenshot(demo_dir, tmp_path):
    QSettings.setDefaultFormat(QSettings.IniFormat)
    QSettings.setPath(QSettings.IniFormat, QSettings.UserScope, str(tmp_path / "settings"))
    app = QApplication.instance() or create_application(["sanger-seek-test"])
    window = MainWindow()
    window.show()
    window.add_paths([str(demo_dir)])

    _wait_until(app, lambda: window._pending == 0)
    assert len(window.project.samples) == 3
    assert window.project.reference is not None
    assert all(sample.analyzed for sample in window.project.samples)

    window.sample_list.setCurrentRow(0)
    app.processEvents()
    assert window.current_sample is not None
    assert window.variants.model.rowCount() >= 3
    first = window.current_sample.variants[0]
    window._select_variant(first)
    assert window.cursor_ref == first.ref_pos
    assert window.qc.table.rowCount() == len(window.current_sample.reads)

    before = window.forward_trace.plot.getViewBox().viewRange()[0]
    window.zoom_traces(0.67)
    app.processEvents()
    after = window.forward_trace.plot.getViewBox().viewRange()[0]
    assert after[1] - after[0] < before[1] - before[0]

    control_files = [
        str(demo_dir / "Sample002_F.ab1"),
        str(demo_dir / "Sample002_F.seq"),
    ]
    window.add_control_paths(control_files)
    _wait_until(app, lambda: window._pending == 0)
    assert window.project.wt_control is not None
    assert len(window.project.wt_control.reads) == 1
    window._select_sample("Sample001")
    app.processEvents()
    visible_traces = [trace for trace in window.trace_views if trace.read is not None]
    assert len(visible_traces) == 3
    sample1 = window.project.sample_by_key("Sample001")
    assert sample1 is not None
    assert any(v.control_status == "absent" for v in sample1.variants)

    # A user pan/zoom on one raw-trace x-axis is propagated through reference
    # coordinates, not raw sample indices, to every other visible trace.
    source = visible_traces[0]
    ref_lo, ref_hi = window.cursor_ref - 12, window.cursor_ref + 12
    source.plot.getViewBox().setXRange(source._ref2x[ref_lo], source._ref2x[ref_hi], padding=0)
    app.processEvents()
    synced = [trace.visible_reference_range() for trace in visible_traces]
    assert all(rng is not None for rng in synced)
    assert max(rng[0] for rng in synced) - min(rng[0] for rng in synced) < 1.5
    assert max(rng[1] for rng in synced) - min(rng[1] for rng in synced) < 1.5

    shot = tmp_path / "desktop-smoke.png"
    window.request_screenshot(str(shot))
    _wait_until(app, shot.exists)
    assert shot.stat().st_size > 10_000

    window._dirty = False
    window.close()


def test_macos_packaging_inputs_exist():
    root = Path(__file__).resolve().parents[1]
    assert (root / "packaging" / "SangerSeek.spec").is_file()
    assert (root / "packaging" / "launcher.py").is_file()
    assert (root / "packaging" / "entitlements.plist").is_file()
    assert (root / "sanger_seek" / "resources" / "app-icon.png").stat().st_size > 10_000
    assert os.access(root / "scripts" / "build_macos.sh", os.X_OK)
