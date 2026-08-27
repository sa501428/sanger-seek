"""Sequence alignment strip: reference, protein and read rows around the
cursor, with mismatch/mixed/indel highlighting. Always centered on the
current cursor position; clicking a column moves the cursor everywhere."""

from __future__ import annotations

from PySide6.QtCore import QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import QSizePolicy, QWidget

from ..core.dna import SETS_TO_IUPAC, bases_match, expand, translate_codon
from ..core.model import Read, Reference, Sample
from .theme import (
    CURSOR_COLOR,
    GAP_BG,
    MISMATCH_BG,
    MIXED_BG,
    base_color,
)

CELL_W = 16
ROW_H = 19
RULER_H = 15
LABEL_W = 116


class AlignmentStrip(QWidget):
    cursorRequested = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.reference: Reference | None = None
        self.sample: Sample | None = None
        self.track_rows: list[tuple[str, Read]] = []
        self.wt_reads: list[Read] = []
        self.sample_reads: list[Read] = []
        self.show_difference = False
        self.cursor: int = 0
        self._ins_map: dict[str, dict[int, str]] = {}
        self._font = QFont("Menlo, Monaco, Courier New", 11)
        self._font_small = QFont("Menlo, Monaco, Courier New", 9)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setMinimumHeight(RULER_H + 3 * ROW_H)
        self.setFocusPolicy(Qt.ClickFocus)

    # ------------------------------------------------------------------ data

    def set_data(
        self,
        reference: Reference | None,
        sample: Sample | None,
        track_rows: list[tuple[str, Read]] | None = None,
        wt_reads: list[Read] | None = None,
        sample_reads: list[Read] | None = None,
    ) -> None:
        self.reference = reference
        self.sample = sample
        self.track_rows = list(track_rows or [])
        self.wt_reads = list(wt_reads or [])
        self.sample_reads = list(sample_reads or [])
        self.show_difference = bool(self.wt_reads and self.sample_reads)
        self._ins_map.clear()
        for _label, read in self.track_rows:
            self._ins_map[read.id] = self._insertions(read)
        rows = 1 + len(self.track_rows) + int(self.show_difference)
        rows += int(reference is not None and bool(reference.cds))
        self.setFixedHeight(RULER_H + max(rows, 2) * ROW_H + 6)
        self.update()

    @staticmethod
    def _insertions(read: Read) -> dict[int, str]:
        """ref junction position -> inserted bases (read insertions)."""
        aln = read.alignment
        if aln is None:
            return {}
        out: dict[int, str] = {}
        rc, fc = aln.read_start, aln.ref_start
        for i, (op, k) in enumerate(aln.ops):
            if op in ("=", "X"):
                rc += k
                fc += k
            elif op == "I":
                if 0 < i < len(aln.ops) - 1:
                    out[fc] = "".join(read.oriented_base(rc + t) for t in range(k))
                rc += k
            else:
                fc += k
        return out

    def set_cursor(self, refpos: int) -> None:
        self.cursor = refpos
        self.update()

    # ---------------------------------------------------------------- paint

    def paintEvent(self, event) -> None:
        p = QPainter(self)
        p.fillRect(self.rect(), QColor("#ffffff"))
        if self.reference is None:
            p.setPen(QColor("#868e96"))
            p.drawText(self.rect(), Qt.AlignCenter, "Load a reference (FASTA / GenBank) to see the alignment")
            p.end()
            return

        ref = self.reference
        n_vis = max((self.width() - LABEL_W) // CELL_W, 5)
        start = self.cursor - n_vis // 2
        fm = QFontMetrics(self._font)

        y = 0
        # ruler
        p.setFont(self._font_small)
        p.setPen(QColor("#868e96"))
        step = max(1, (n_vis // 6 // 5) * 5 or 5)
        for col in range(n_vis):
            g = start + col
            if 0 <= g < ref.n and (g + 1) % step == 0:
                x = LABEL_W + col * CELL_W
                p.drawText(QRectF(x - 30, y, CELL_W + 60, RULER_H), Qt.AlignCenter, f"{g + 1}")
        y += RULER_H

        has_cds = bool(ref.cds)
        if has_cds:
            self._paint_protein_row(p, y, start, n_vis)
            y += ROW_H

        self._paint_row_label(p, y, "Reference")
        p.setFont(self._font)
        for col in range(n_vis):
            g = start + col
            if not (0 <= g < ref.n):
                continue
            base = ref.seq[g]
            self._cell(p, col, y, base, base_color(base), None)
        y += ROW_H

        for label, read in self.track_rows:
            self._paint_read_row(p, y, read, start, n_vis, label)
            y += ROW_H

        if self.show_difference:
            self._paint_difference_row(p, y, start, n_vis)
            y += ROW_H

        # cursor column outline
        col = self.cursor - start
        if 0 <= col < n_vis:
            p.setPen(QPen(QColor(CURSOR_COLOR), 1.6))
            p.setBrush(Qt.NoBrush)
            p.drawRect(QRectF(LABEL_W + col * CELL_W, RULER_H, CELL_W, y - RULER_H))
        p.end()
        del fm

    def _paint_row_label(self, p: QPainter, y: int, text: str, color: str = "#495057") -> None:
        p.setFont(self._font_small)
        p.setPen(QColor(color))
        p.drawText(QRectF(4, y, LABEL_W - 8, ROW_H), Qt.AlignVCenter | Qt.AlignLeft, text)

    def _cell(self, p: QPainter, col: int, y: int, text: str, fg: str, bg: str | None) -> None:
        x = LABEL_W + col * CELL_W
        r = QRectF(x, y, CELL_W, ROW_H)
        if bg:
            p.fillRect(r, QColor(bg))
        p.setFont(self._font)
        p.setPen(QColor(fg))
        p.drawText(r, Qt.AlignCenter, text)

    def _paint_protein_row(self, p: QPainter, y: int, start: int, n_vis: int) -> None:
        ref = self.reference
        self._paint_row_label(p, y, "Protein")
        feature = None
        for col in range(n_vis):
            g = start + col
            if not (0 <= g < ref.n):
                continue
            feature = ref.feature_at(g)
            if feature is None:
                continue
            ci = feature.cds_index(g)
            if ci is None:
                continue
            cpos = ci - (feature.codon_start - 1)
            if cpos < 0:
                continue
            codon_i, within = divmod(cpos, 3)
            # light alternating codon shading
            if codon_i % 2 == 0:
                p.fillRect(QRectF(LABEL_W + col * CELL_W, y, CELL_W, ROW_H), QColor("#f1f3f5"))
            if within == 1:  # middle base: draw the aa letter
                coding = feature.cds_seq[feature.codon_start - 1 :]
                codon = coding[codon_i * 3 : codon_i * 3 + 3]
                if len(codon) == 3:
                    aa = translate_codon(codon)
                    p.setFont(self._font)
                    p.setPen(QColor("#495057"))
                    p.drawText(
                        QRectF(LABEL_W + col * CELL_W, y, CELL_W, ROW_H),
                        Qt.AlignCenter,
                        aa,
                    )

    def _paint_read_row(
        self, p: QPainter, y: int, read: Read, start: int, n_vis: int, label: str | None = None
    ) -> None:
        label = label or (
            "Forward" if read.orientation == "F" else "Reverse" if read.orientation == "R" else read.label[:14]
        )
        self._paint_row_label(p, y, label)
        aln = read.alignment
        if aln is None:
            p.setFont(self._font_small)
            p.setPen(QColor("#adb5bd"))
            p.drawText(
                QRectF(LABEL_W, y, self.width() - LABEL_W, ROW_H),
                Qt.AlignVCenter | Qt.AlignLeft,
                read.error or "not aligned",
            )
            return
        ref = self.reference
        ins_map = self._ins_map.get(read.id, {})
        for col in range(n_vis):
            g = start + col
            if not (0 <= g < ref.n) or not (aln.ref_start <= g < aln.ref_end):
                continue
            idx = aln.read_index_at(g)
            if idx is None:
                self._cell(p, col, y, "–", "#495057", GAP_BG)
            else:
                call = read.oriented_base(idx)
                ref_base = ref.seq[g]
                if call == ref_base:
                    bg = None
                elif ref_base in expand(call) and len(expand(call)) > 1:
                    bg = MIXED_BG
                else:
                    bg = MISMATCH_BG
                self._cell(p, col, y, call, base_color(call), bg)
            if g in ins_map:
                x = LABEL_W + col * CELL_W
                p.setPen(QPen(QColor("#9c36b5"), 2))
                p.drawLine(int(x), int(y + 2), int(x), int(y + ROW_H - 2))

    @staticmethod
    def _group_call(reads: list[Read], refpos: int) -> str | None:
        """Return a compact IUPAC consensus for the reads covering refpos."""
        calls: list[str] = []
        for read in reads:
            aln = read.alignment
            if aln is None or not (aln.ref_start <= refpos < aln.ref_end):
                continue
            idx = aln.read_index_at(refpos)
            calls.append("-" if idx is None else read.oriented_base(idx))
        if not calls:
            return None
        if all(call == "-" for call in calls):
            return "-"
        concrete: set[str] = set()
        for call in calls:
            concrete.update(expand(call))
        return SETS_TO_IUPAC.get(frozenset(concrete), "N") if concrete else "-"

    def difference_call_at(self, refpos: int) -> tuple[str | None, str | None]:
        """WT and assessed consensus calls used by the subtraction row."""
        return self._group_call(self.wt_reads, refpos), self._group_call(self.sample_reads, refpos)

    def _paint_difference_row(self, p: QPainter, y: int, start: int, n_vis: int) -> None:
        self._paint_row_label(p, y, "Sample − WT", "#9c36b5")
        for col in range(n_vis):
            g = start + col
            wt, sample = self.difference_call_at(g)
            if wt is None or sample is None:
                continue
            if wt == sample:
                self._cell(p, col, y, "·", "#ced4da", None)
                continue
            overlaps = wt != "-" and sample != "-" and bases_match(wt, sample)
            bg = MIXED_BG if overlaps else MISMATCH_BG
            self._cell(p, col, y, sample, base_color(sample), bg)

    # ---------------------------------------------------------------- input

    def mousePressEvent(self, event) -> None:
        if self.reference is None:
            return
        n_vis = max((self.width() - LABEL_W) // CELL_W, 5)
        start = self.cursor - n_vis // 2
        col = int((event.position().x() - LABEL_W) // CELL_W)
        if col >= 0:
            g = start + col
            if 0 <= g < self.reference.n:
                self.cursorRequested.emit(g)

    def wheelEvent(self, event) -> None:
        if self.reference is None:
            return
        delta = event.angleDelta().y() or event.angleDelta().x()
        step = -3 if delta > 0 else 3
        g = max(0, min(self.reference.n - 1, self.cursor + step))
        if g != self.cursor:
            self.cursorRequested.emit(g)
