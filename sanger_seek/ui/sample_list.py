"""Cases dock: one row per case with read/variant/QC summary."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QListWidget, QListWidgetItem

from ..core.model import Sample


class SampleList(QListWidget):
    sampleSelected = Signal(str)     # sample key

    def __init__(self, parent=None):
        super().__init__(parent)
        self.currentItemChanged.connect(self._on_current)

    def set_samples(self, samples: list[Sample], keep_key: str | None = None) -> None:
        self.blockSignals(True)
        self.clear()
        current_row = 0
        for i, s in enumerate(samples):
            item = QListWidgetItem(self._text(s))
            item.setData(Qt.UserRole, s.key)
            if s.error:
                item.setToolTip(s.error)
            self.addItem(item)
            if s.key == keep_key:
                current_row = i
        self.blockSignals(False)
        if self.count():
            self.setCurrentRow(current_row)

    @staticmethod
    def _text(s: Sample) -> str:
        n_f = sum(1 for r in s.reads if r.orientation == "F")
        n_r = sum(1 for r in s.reads if r.orientation == "R")
        reads = f"{len(s.reads)} read{'s' if len(s.reads) != 1 else ''}"
        if n_f or n_r:
            reads += f" ({n_f}F/{n_r}R)"
        if not s.analyzed:
            status = "analyzing…" if not s.error else f"error: {s.error}"
            return f"{s.name}\n{reads} · {status}"
        n_high = sum(1 for v in s.variants if v.confidence == "high")
        n_rev = sum(1 for v in s.variants if v.confidence == "review")
        vs = f"{len(s.variants)} variant{'s' if len(s.variants) != 1 else ''}"
        if s.variants:
            vs += f" ({n_high} high, {n_rev} review)"
        warn = " ⚠ review" if s.qc_flags or any(r.error or r.qc_flags for r in s.reads) else ""
        return f"{s.name}{warn}\n{reads} · {vs}"

    def refresh_sample(self, sample: Sample) -> None:
        for i in range(self.count()):
            item = self.item(i)
            if item.data(Qt.UserRole) == sample.key:
                item.setText(self._text(sample))
                return

    def _on_current(self, current, previous) -> None:
        if current is not None:
            self.sampleSelected.emit(current.data(Qt.UserRole))
