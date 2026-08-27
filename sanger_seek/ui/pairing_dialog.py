"""Manual case grouping and forward/reverse assignment."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from ..core.model import Read, Sample


class PairingDialog(QDialog):
    """Edit case membership and orientation hints for all reads."""

    def __init__(self, samples: list[Sample], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Assign Reads to Cases")
        self.resize(720, 420)
        self._rows: list[Read] = []

        layout = QVBoxLayout(self)
        help_text = QLabel(
            "Give matching forward and reverse reads the same case name. "
            "Choose a role when filenames are ambiguous; Auto lets alignment decide."
        )
        help_text.setWordWrap(True)
        layout.addWidget(help_text)

        reads = [(sample, read) for sample in samples for read in sample.reads]
        self.table = QTableWidget(len(reads), 4)
        self.table.setHorizontalHeaderLabels(["Read", "Current files", "Case", "Role"])
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)

        for row, (sample, read) in enumerate(reads):
            self._rows.append(read)
            label = QTableWidgetItem(read.label)
            label.setFlags(label.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(row, 0, label)
            files = " / ".join(
                name for name in (
                    read.ab1_path and read.ab1_path.rsplit("/", 1)[-1],
                    read.seq_path and read.seq_path.rsplit("/", 1)[-1],
                    read.phd_path and read.phd_path.rsplit("/", 1)[-1],
                ) if name
            )
            file_item = QTableWidgetItem(files)
            file_item.setFlags(file_item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(row, 1, file_item)
            sample_edit = QLineEdit(sample.name)
            sample_edit.setProperty("assignmentSample", True)
            self.table.setCellWidget(row, 2, sample_edit)
            role = QComboBox()
            role.addItem("Auto", None)
            role.addItem("Forward", "F")
            role.addItem("Reverse", "R")
            current = read.orientation_override or read.orientation_hint or (
                read.orientation if read.orientation in ("F", "R") else None
            )
            role.setCurrentIndex({None: 0, "F": 1, "R": 2}[current])
            self.table.setCellWidget(row, 3, role)
        layout.addWidget(self.table, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Cancel | QDialogButtonBox.Save)
        buttons.rejected.connect(self.reject)
        buttons.accepted.connect(self._validate_and_accept)
        layout.addWidget(buttons)

    def _validate_and_accept(self) -> None:
        if any(not self.table.cellWidget(row, 2).text().strip() for row in range(self.table.rowCount())):
            QMessageBox.warning(self, "Case name required", "Every read must have a case name.")
            return
        self.accept()

    def assignments(self) -> list[tuple[Read, str, str | None]]:
        out = []
        for row, read in enumerate(self._rows):
            sample_name = self.table.cellWidget(row, 2).text().strip()
            orientation = self.table.cellWidget(row, 3).currentData()
            out.append((read, sample_name, orientation))
        return out
