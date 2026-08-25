"""Top summary strip: sample, reference, reads, variants, QC toggle."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QWidget

from ..core.model import Reference, Sample


def _chip(text: str = "") -> QLabel:
    lbl = QLabel(text)
    lbl.setProperty("chip", "true")
    return lbl


class SummaryBar(QWidget):
    qcToggled = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setSpacing(8)

        self.sample_label = QLabel("No sample")
        self.sample_label.setStyleSheet("font-weight: 700; font-size: 14px;")
        layout.addWidget(self.sample_label)

        self.ref_chip = _chip("no reference")
        self.reads_chip = _chip("0 reads")
        self.var_chip = _chip("–")
        self.control_chip = _chip("WT: not loaded")
        self.disc_chip = _chip("")
        for c in (self.ref_chip, self.control_chip, self.reads_chip, self.var_chip, self.disc_chip):
            layout.addWidget(c)
        self.disc_chip.hide()

        layout.addStretch(1)
        self.qc_btn = QPushButton("QC details")
        self.qc_btn.setCheckable(True)
        self.qc_btn.toggled.connect(self.qcToggled.emit)
        layout.addWidget(self.qc_btn)

    def update_summary(
        self, sample: Sample | None, reference: Reference | None, control: Sample | None = None
    ) -> None:
        if reference is not None:
            gene = reference.cds[0].gene if reference.cds else None
            self.ref_chip.setText(
                f"ref: {reference.name}" + (f" · {gene}" if gene else "") + f" · {reference.n} bp"
            )
        else:
            self.ref_chip.setText("no reference — load FASTA/GenBank")

        if control is None:
            self.control_chip.setText("WT: not loaded")
        elif control.analyzed:
            self.control_chip.setText(f"WT: {len(control.reads)} contig read(s)")
        else:
            self.control_chip.setText("WT: analyzing…")

        if sample is None:
            self.sample_label.setText("No sample")
            self.reads_chip.setText("0 reads")
            self.var_chip.setText("–")
            self.disc_chip.hide()
            return

        self.sample_label.setText(sample.name)
        n_f = len(sample.forward_reads)
        n_r = len(sample.reverse_reads)
        self.reads_chip.setText(f"{len(sample.reads)} reads ({n_f}F/{n_r}R)")
        if not sample.analyzed:
            self.var_chip.setText("analyzing…")
        else:
            n_high = sum(1 for v in sample.variants if v.confidence == "high")
            self.var_chip.setText(
                f"{len(sample.variants)} candidate variants · {n_high} high confidence"
            )
        n_disc = sum(r.discrepancies.count for r in sample.reads if r.discrepancies)
        if n_disc:
            self.disc_chip.setText(f"⚠ {n_disc} .seq/.ab1 call difference(s)")
            self.disc_chip.show()
        else:
            self.disc_chip.hide()
