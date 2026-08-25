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
