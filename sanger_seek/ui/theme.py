"""Colors and application-wide styling."""

from __future__ import annotations

# Classic chromatogram channel colors
BASE_COLORS = {
    "A": "#2f9e44",   # green
    "C": "#1971c2",   # blue
    "G": "#212529",   # near-black
    "T": "#e03131",   # red
}
AMBIG_COLOR = "#9c36b5"

CONFIDENCE_COLORS = {
    "high": "#2b8a3e",
    "review": "#e8890c",
    "low": "#868e96",
}
CONFIDENCE_BG = {
    "high": "#ebfbee",
    "review": "#fff4e6",
    "low": "#f1f3f5",
}

MISMATCH_BG = "#ffe3e3"
MIXED_BG = "#fff3bf"
GAP_BG = "#e9ecef"
CURSOR_COLOR = "#7048e8"
QUAL_COLORS = {"good": "#a5d8ff", "mid": "#ffd8a8", "bad": "#ffc9c9"}

APP_QSS = """
QMainWindow, QDialog { background: #f8f9fa; }
QWidget { font-size: 12px; color: #212529; }
QToolBar { background: #ffffff; border-bottom: 1px solid #dee2e6; spacing: 6px; padding: 3px; }
QStatusBar { background: #ffffff; border-top: 1px solid #dee2e6; }
QDockWidget::title { background: #f1f3f5; padding: 5px; font-weight: bold; }
QListWidget, QTableView, QTextBrowser {
    background: #ffffff; border: 1px solid #dee2e6; border-radius: 4px;
}
QTableView { gridline-color: #f1f3f5; selection-background-color: #e7f5ff; selection-color: #212529; }
QHeaderView::section {
    background: #f8f9fa; border: none; border-bottom: 1px solid #dee2e6;
    border-right: 1px solid #f1f3f5; padding: 4px 6px; font-weight: 600;
}
QPushButton, QComboBox {
    background: #ffffff; border: 1px solid #ced4da; border-radius: 4px; padding: 4px 10px;
}
QPushButton:hover, QComboBox:hover { border-color: #74c0fc; }
QPushButton:checked { background: #e7f5ff; border-color: #4dabf7; }
QSplitter::handle { background: #dee2e6; }
QSplitter::handle:vertical { height: 3px; }
QSplitter::handle:horizontal { width: 3px; }
QLabel[chip="true"] {
    background: #f1f3f5; border: 1px solid #dee2e6; border-radius: 8px; padding: 2px 8px;
}
"""


def base_color(base: str) -> str:
    return BASE_COLORS.get(base.upper(), AMBIG_COLOR)
