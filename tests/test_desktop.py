"""Offscreen smoke coverage for the assembled desktop application."""

import os
import time
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QPointF, QSettings, Qt  # noqa: E402
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
        str(demo_dir / "Sample002_R.ab1"),
        str(demo_dir / "Sample002_R.seq"),
    ]
    window.add_control_paths(control_files)
    _wait_until(app, lambda: window._pending == 0)
    assert window.project.wt_control is not None
    assert len(window.project.wt_control.reads) == 2
    window._select_sample("Sample001")
    app.processEvents()
    visible_traces = [trace for trace in window.trace_views if trace.read is not None]
    assert len(visible_traces) == 4
    assert [label for label, _read in window.alignment.track_rows] == [
        "WT Fwd", "WT Rev", "Sample Fwd", "Sample Rev",
    ]
    assert window.alignment.show_difference
    assert len(window.alignment.track_rows) + 1 == 5  # reference + one sequence per trace
    wt_call, sample_call = window.alignment.difference_call_at(window.cursor_ref)
    assert wt_call is not None and sample_call is not None
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
    assert window.cursor_ref == round((ref_lo + ref_hi) / 2)

    # Double-clicking a base selects it and requests one 2x synchronized zoom.
    clicked_ref = window.cursor_ref + 2
    clicked_scene_pos = source.plot.getViewBox().mapViewToScene(
        QPointF(source._ref2x[clicked_ref], 0)
    )
    before_span = source.visible_reference_range()[1] - source.visible_reference_range()[0]
    event = SimpleNamespace(
        button=lambda: Qt.LeftButton,
        scenePos=lambda: clicked_scene_pos,
        double=lambda: True,
    )
    source._on_click(event)
    app.processEvents()
    assert window.cursor_ref == clicked_ref
    zoomed = [trace.visible_reference_range() for trace in visible_traces]
    assert all(rng is not None for rng in zoomed)
    assert source.visible_reference_range()[1] - source.visible_reference_range()[0] < before_span * 0.6
    assert all(abs((rng[0] + rng[1]) / 2 - clicked_ref) < 1 for rng in zoomed)

    assignments = []
    for sample in window.project.samples:
        for read in sample.reads:
            name = "ManualPair" if sample.key == "Sample001" else sample.name
            role = "F" if read.label.endswith("_F") else "R" if read.label.endswith("_R") else None
            assignments.append((read, name, role))
    window.apply_read_assignments(assignments)
    _wait_until(app, lambda: window._pending == 0)
    manual = window.project.sample_by_key("ManualPair")
    assert manual is not None and len(manual.reads) == 2
    assert {read.orientation_override for read in manual.reads} == {"F", "R"}

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
