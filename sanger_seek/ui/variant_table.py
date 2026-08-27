"""Variant review table with visualization filters."""

from __future__ import annotations

from PySide6.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    QSortFilterProxyModel,
    Qt,
    Signal,
)
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QTableView,
    QVBoxLayout,
    QWidget,
)

from ..core.model import Variant
from .theme import CONFIDENCE_BG, CONFIDENCE_COLORS

COLUMNS = ["Variant", "Protein", "Effect", "Fwd", "Rev", "Trace", "Confidence", "Mix"]

FILTERS = [
    ("All differences", lambda v: True),
    ("Coding only", lambda v: v.effect not in (None, "non-coding")),
    ("Missense / nonsense", lambda v: v.effect in ("missense", "nonsense", "stop-loss")),
    ("Indels", lambda v: v.kind in ("ins", "del")),
    ("Mixed peaks", lambda v: v.mixed),
    ("Strand disagreement", lambda v: v.strand_status == "discordant"),
    ("Needs review", lambda v: v.confidence in ("review", "low")),
]


def strand_mark(v: Variant, orientation: str) -> str:
    evs = [e for e in v.evidence if e.orientation == orientation]
    covered = [e for e in evs if e.covered]
    if not covered:
        return "–"
    if any(e.supports for e in covered):
        return "✓"
    if all(e.supports is None for e in covered):
        return "?"
    return "✗"


class VariantModel(QAbstractTableModel):
    def __init__(self):
        super().__init__()
        self.variants: list[Variant] = []

    def set_variants(self, variants: list[Variant]) -> None:
        self.beginResetModel()
        self.variants = variants
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.variants)

    def columnCount(self, parent=QModelIndex()):
        return len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return COLUMNS[section]
        return None

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        v = self.variants[index.row()]
        col = index.column()
        if role == Qt.DisplayRole:
            if col == 0:
                return v.label
            if col == 1:
                return v.protein or ""
            if col == 2:
                return v.effect or ""
            if col == 3:
                return strand_mark(v, "F")
            if col == 4:
                return strand_mark(v, "R")
            if col == 5:
                return v.trace_quality.capitalize()
            if col == 6:
                return v.confidence.capitalize()
            if col == 7:
                if v.mixed and v.alt_fraction is not None:
                    ref_disp = v.ref_bases or "·"
                    return f"{ref_disp} {100 - v.alt_fraction * 100:.0f}% / {v.alt_bases or '·'} {v.alt_fraction * 100:.0f}%"
                return "mixed" if v.mixed else ""
        if role == Qt.TextAlignmentRole and col in (3, 4):
            return Qt.AlignCenter
        if role == Qt.BackgroundRole:
            return QBrush(QColor(CONFIDENCE_BG.get(v.confidence, "#ffffff")))
        if role == Qt.ForegroundRole and col == 6:
            return QBrush(QColor(CONFIDENCE_COLORS.get(v.confidence, "#212529")))
        if role == Qt.ToolTipRole:
            return self._tooltip(v)
        return None

    @staticmethod
    def _tooltip(v: Variant) -> str:
        lines = [v.g_label if v.cdna else v.label]
        if v.gene:
            lines.append(f"gene: {v.gene}")
        lines.append(f"strand status: {v.strand_status}")
        for e in v.evidence:
            if not e.covered:
                lines.append(f"{e.read_label}: not covered")
                continue
            s = "supports" if e.supports else "no support" if e.supports is False else "ambiguous"
            q = f" Q{e.qual}" if e.qual is not None else ""
            r = f" 2nd/1st {e.ratio:.2f}" if e.ratio is not None else ""
            lines.append(f"{e.read_label}: {e.call or '?'}{q}{r} — {s}")
        return "\n".join(lines)

    def variant_at(self, row: int) -> Variant | None:
        return self.variants[row] if 0 <= row < len(self.variants) else None


class VariantFilterProxy(QSortFilterProxyModel):
    def __init__(self):
        super().__init__()
        self.mode = 0

    def set_mode(self, mode: int) -> None:
        self.mode = mode
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        v = self.sourceModel().variant_at(source_row)
        if v is None:
            return False
        return FILTERS[self.mode][1](v)


class VariantPanel(QWidget):
    variantSelected = Signal(object)     # Variant

    def __init__(self, parent=None):
        super().__init__(parent)
        self.model = VariantModel()
        self.proxy = VariantFilterProxy()
        self.proxy.setSourceModel(self.model)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)

        top = QHBoxLayout()
        title = QLabel("Variants")
        title.setStyleSheet("font-weight: 600; color: #495057;")
        top.addWidget(title)
        self.filter_box = QComboBox()
        self.filter_box.addItems([name for name, _ in FILTERS])
        self.filter_box.currentIndexChanged.connect(self._on_filter)
        top.addWidget(self.filter_box)
        self.count_label = QLabel("")
        self.count_label.setMinimumWidth(42)
        self.count_label.setStyleSheet("color: #868e96;")
        top.addWidget(self.count_label)
        top.addStretch(1)
        layout.addLayout(top)

        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSelectionBehavior(QTableView.SelectRows)
        self.table.setSelectionMode(QTableView.SingleSelection)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setAlternatingRowColors(False)
        self.table.clicked.connect(self._on_click)
        layout.addWidget(self.table, 1)

    def set_variants(self, variants: list[Variant]) -> None:
        self.model.set_variants(variants)
        self._update_count()

    def _on_filter(self, mode: int) -> None:
        self.proxy.set_mode(mode)
        self._update_count()

    def _update_count(self) -> None:
        shown = self.proxy.rowCount()
        total = self.model.rowCount()
        self.count_label.setText(f"{shown} of {total}" if shown != total else f"{total}")

    def _on_click(self, proxy_index) -> None:
        src = self.proxy.mapToSource(proxy_index)
        v = self.model.variant_at(src.row())
        if v is not None:
            self.variantSelected.emit(v)

    def select_variant(self, variant: Variant) -> None:
        for row in range(self.model.rowCount()):
            if self.model.variant_at(row) is variant:
                src_idx = self.model.index(row, 0)
                p_idx = self.proxy.mapFromSource(src_idx)
                if p_idx.isValid():
                    self.table.selectRow(p_idx.row())
                return

    def visible_variants(self) -> list[Variant]:
        out = []
        for prow in range(self.proxy.rowCount()):
            src = self.proxy.mapToSource(self.proxy.index(prow, 0))
            v = self.model.variant_at(src.row())
            if v is not None:
                out.append(v)
        return out
