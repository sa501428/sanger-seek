"""Progressively disclosed read and per-position quality details."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFormLayout,
    QGroupBox,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..core.model import Read, Sample


class QCPanel(QWidget):
    """Shows auditable QC values for the current reference position."""

    HEADERS = ["Read", "Call", "Q", "Primary", "Secondary", "Ratio", "Noise", "Trim edge"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.sample: Sample | None = None
        self.refpos = 0

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        self.position_label = QLabel("Select a reference position")
        self.position_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.position_label)

        self.table = QTableWidget(0, len(self.HEADERS))
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.table, 1)

        imported = QGroupBox("Imported call comparison")
        form = QFormLayout(imported)
        self.discrepancy = QLabel("No sample")
        self.discrepancy.setWordWrap(True)
        form.addRow(self.discrepancy)
        layout.addWidget(imported)

    def set_data(self, sample: Sample | None, refpos: int) -> None:
        self.sample = sample
        self.refpos = refpos
        self.position_label.setText(f"Reference position {refpos + 1}")
        self.table.setRowCount(0)
        if sample is None:
            self.discrepancy.setText("No sample")
            return

        for read in sample.reads:
            self._add_read(read, refpos)
        reports = []
        for read in sample.reads:
            report = read.discrepancies
            if report is None:
                continue
            if report.count == 0:
                reports.append(f"{read.label}: AB1 and .seq calls agree")
            else:
                where = ", ".join(str(p + 1) for p in report.positions[:12])
                if len(report.positions) > 12:
                    where += ", …"
                detail = f" at AB1 call(s) {where}" if where else ""
                if report.note:
                    detail += f"; {report.note}"
                reports.append(f"{read.label}: {report.count} difference(s){detail}")
        self.discrepancy.setText("\n".join(reports) if reports else "No paired AB1/.seq calls")

    def _add_read(self, read: Read, refpos: int) -> None:
        aln = read.alignment
        idx = aln.read_index_at(refpos) if aln is not None else None
        values: list[str]
        if idx is None:
            values = [read.label, "—", "—", "—", "—", "—", "—", "—"]
        else:
            peak = read.peak_at(idx)
            values = [
                read.label,
                read.oriented_base(idx),
                str(read.qual_at(idx)) if read.qual_at(idx) is not None else "—",
                f"{peak.primary_h:.0f}" if peak else "—",
                (f"{peak.secondary_base} {peak.secondary_h:.0f}" if peak else "—"),
                f"{peak.ratio:.2f}" if peak else "—",
                f"{peak.noise:.0f}" if peak else "—",
                str(read.dist_from_trim_end(idx)),
            ]
        row = self.table.rowCount()
        self.table.insertRow(row)
        for col, value in enumerate(values):
            item = QTableWidgetItem(value)
            if col > 0:
                item.setTextAlignment(Qt.AlignCenter)
            self.table.setItem(row, col, item)
