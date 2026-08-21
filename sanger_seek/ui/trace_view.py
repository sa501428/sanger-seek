"""Interactive chromatogram view (pyqtgraph).

Shows one read's raw traces in *reference orientation* (reverse reads are
mirrored with complemented channels), with base calls above peaks, quality
bars, a reference-coordinate ruler, variant flags, and a cursor crosshair.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

from ..core.dna import complement, revcomp
from ..core.model import Read, Reference, Variant
from .theme import CONFIDENCE_COLORS, CURSOR_COLOR, base_color

pg.setConfigOptions(antialias=True, background="w", foreground="#343a40")

MAX_BASE_LABELS = 160


class RefAxis(pg.AxisItem):
    """Bottom axis labeled in 1-based reference coordinates."""

    def __init__(self):
        super().__init__(orientation="bottom")
        self._xs = np.zeros(0)          # peak x positions (ascending)
        self._refs = np.zeros(0)        # ref pos per peak (-1 = unaligned)
        self._ref2x: dict[int, float] = {}

    def set_mapping(self, xs: np.ndarray, refs: np.ndarray) -> None:
        self._xs = xs
        self._refs = refs
        self._ref2x = {int(r): float(x) for x, r in zip(xs, refs) if r >= 0}
        self.picture = None
        self.update()

    def tickValues(self, minVal, maxVal, size):
        if not self._ref2x:
            return []
        i0, i1 = np.searchsorted(self._xs, [minVal, maxVal])
        visible = self._refs[max(i0 - 1, 0) : i1 + 1]
        visible = visible[visible >= 0]
        if len(visible) < 2:
            return []
        rmin, rmax = int(visible.min()), int(visible.max())
        span = max(rmax - rmin, 1)
        target = max(int(size / 90), 1)
        raw = span / target
        mag = 10 ** int(np.floor(np.log10(raw))) if raw > 0 else 1
        step = next((s * mag for s in (1, 2, 5, 10) if s * mag >= raw), 10 * mag)
        step = max(int(step), 1)
        first = ((rmin + step) // step) * step
        major = [
            self._ref2x[r]
            for r in range(first, rmax + 1, step)
            if r in self._ref2x
        ]
        return [(step, major)]

    def tickStrings(self, values, scale, spacing):
        out = []
        for v in values:
            i = int(np.searchsorted(self._xs, v))
            i = min(max(i, 0), len(self._refs) - 1)
            # snap to the peak nearest this x
            if i > 0 and abs(self._xs[i - 1] - v) < abs(self._xs[i] - v):
                i -= 1
            r = int(self._refs[i])
            out.append(str(r + 1) if r >= 0 else "")
        return out


class TraceView(QWidget):
    positionClicked = Signal(int, object)   # refpos (or -1), read
    rangeChanged = Signal()

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self._title = title
        self.read: Read | None = None
        self.reference: Reference | None = None
        self._xs = np.zeros(0)
        self._refs = np.zeros(0)
        self._calls = ""
        self._quals: np.ndarray | None = None
        self._ref2x: dict[int, float] = {}
        self._ymax = 1.0
        self._labels: list[pg.TextItem] = []
        self._marker_items: list = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        self.header = QLabel(title)
        self.header.setStyleSheet("font-weight: 600; color: #495057; padding: 1px 4px;")
        layout.addWidget(self.header)

        self.axis = RefAxis()
        self.plot = pg.PlotWidget(axisItems={"bottom": self.axis})
        self.plot.setMouseEnabled(x=True, y=False)
        self.plot.hideButtons()
        self.plot.setMenuEnabled(False)
        self.plot.getAxis("left").setWidth(10)
        self.plot.getAxis("left").setStyle(showValues=False)
        layout.addWidget(self.plot, 1)

        self.empty = QLabel("No read loaded")
        self.empty.setAlignment(Qt.AlignCenter)
        self.empty.setStyleSheet("color: #868e96; font-style: italic;")
        layout.addWidget(self.empty, 1)
        self.empty.hide()

        self.curves = {}
        for b in "GATC":
            c = pg.PlotCurveItem(pen=pg.mkPen(base_color(b), width=1.4))
            c.setDownsampling(auto=True)
            c.setClipToView(True)
            self.plot.addItem(c)
            self.curves[b] = c

        self.qual_bars: list[pg.BarGraphItem] = []
        self.cursor_line = pg.InfiniteLine(
            angle=90, movable=False,
            pen=pg.mkPen(CURSOR_COLOR, width=1.6, style=Qt.DashLine),
        )
        self.cursor_line.setZValue(50)
        self.plot.addItem(self.cursor_line)
        self.cursor_line.hide()

        self._label_font = QFont("Menlo, Monaco, Courier New", 11)
        self._label_font.setBold(True)

        self._relabel_timer = QTimer(self)
        self._relabel_timer.setSingleShot(True)
        self._relabel_timer.setInterval(25)
        self._relabel_timer.timeout.connect(self._update_base_labels)
        self.plot.getViewBox().sigXRangeChanged.connect(
            lambda *_: (self._relabel_timer.start(), self.rangeChanged.emit())
        )
        self.plot.scene().sigMouseClicked.connect(self._on_click)

    # ------------------------------------------------------------------ data

    def set_read(self, read: Read | None, reference: Reference | None) -> None:
        self.read = read
        self.reference = reference
        for item in self._marker_items:
            self.plot.removeItem(item)
        self._marker_items.clear()
        for bar in self.qual_bars:
            self.plot.removeItem(bar)
        self.qual_bars.clear()

        if read is None or read.trace is None:
            self.plot.hide()
            self.empty.setText(
                "No trace"
                if read is None
                else f"{read.label}: sequence only (no .ab1 chromatogram)"
            )
            self.empty.show()
            self.header.setText(self._title)
            return
        self.plot.show()
        self.empty.hide()

        trace = read.trace
        S = trace.n_samples
        rev = read.is_reverse
        for b in "GATC":
            src = trace.channels.get(complement(b) if rev else b)
            if src is None:
                src = np.zeros(S, dtype=np.int32)
            y = src[::-1] if rev else src
            self.curves[b].setData(np.arange(len(y)), y.astype(np.float64))

        ploc = np.asarray(trace.ploc)
        self._xs = (S - 1 - ploc)[::-1].astype(np.float64) if rev else ploc.astype(np.float64)
        self._calls = revcomp(trace.calls) if rev else trace.calls
        self._quals = trace.quals[::-1] if rev else trace.quals

        n = len(self._calls)
        if read.alignment is not None:
            self._refs = read.alignment.read_to_ref.astype(np.int64)
        else:
            self._refs = np.full(n, -1, dtype=np.int64)
        self.axis.set_mapping(self._xs, self._refs)
        self._ref2x = {int(r): float(x) for x, r in zip(self._xs, self._refs) if r >= 0}

        peak = max(float(np.percentile(
            np.concatenate([c.yData for c in self.curves.values()]), 99.8)), 10.0)
        self._ymax = peak * 1.30
        vb = self.plot.getViewBox()
        vb.setLimits(xMin=0, xMax=S, yMin=0, yMax=self._ymax)
        self.plot.setYRange(0, self._ymax, padding=0)

        self._build_quality_bars()
        self._describe_header()
        self._update_base_labels()

    def _describe_header(self) -> None:
        r = self.read
        parts = [f"{self._title} — {r.label}"]
        if r.quals is not None and len(r.quals):
            parts.append(f"mean Q {float(np.mean(r.quals)):.0f}")
        ts, te = r.trim
        parts.append(f"{r.n} bp (trimmed to {ts + 1}–{te})")
        if r.alignment is not None:
            a = r.alignment
            parts.append(f"ref {a.ref_start + 1}–{a.ref_end} · identity {a.identity * 100:.1f}%")
        else:
            parts.append("not aligned")
        if r.discrepancies is not None and r.discrepancies.count:
            parts.append(f"⚠ {r.discrepancies.count} .seq/.ab1 difference(s)")
        self.header.setText("   ·   ".join(parts))

    def _build_quality_bars(self) -> None:
        if self._quals is None or not len(self._xs):
            return
        q = np.asarray(self._quals, dtype=np.float64)
        y0 = self._ymax * 0.86
        h = (np.clip(q, 0, 60) / 60.0) * self._ymax * 0.105
        spacing = float(np.median(np.diff(self._xs))) if len(self._xs) > 1 else 10.0
        groups = {
            "#74c0fc": q >= 30,
            "#ffc078": (q >= 20) & (q < 30),
            "#ff8787": q < 20,
        }
        for color, mask in groups.items():
            if not mask.any():
                continue
            bar = pg.BarGraphItem(
                x=self._xs[mask], y0=y0, height=h[mask],
                width=spacing * 0.7, brush=pg.mkBrush(color), pen=pg.mkPen(None),
            )
            bar.setOpacity(0.75)
            bar.setZValue(-5)
            self.plot.addItem(bar)
            self.qual_bars.append(bar)

    def set_variants(self, variants: list[Variant]) -> None:
        for item in self._marker_items:
            self.plot.removeItem(item)
        self._marker_items.clear()
        if self.read is None or self.read.alignment is None:
            return
        for v in variants:
            x = self._variant_x(v)
            if x is None:
                continue
            color = CONFIDENCE_COLORS.get(v.confidence, "#868e96")
            line = pg.InfiniteLine(
                pos=x, angle=90, movable=False,
                pen=pg.mkPen(color, width=1.2, style=Qt.DotLine),
            )
            line.setZValue(20)
            label = pg.InfLineLabel(
                line, v.label, position=0.97, color=color, movable=False,
                fill=pg.mkBrush(255, 255, 255, 210),
            )
            label.setZValue(21)
            self.plot.addItem(line)
            self._marker_items.append(line)

    def _variant_x(self, v: Variant) -> float | None:
        if v.kind == "ins":
            xa = self._ref2x.get(v.ref_pos - 1)
            xb = self._ref2x.get(v.ref_pos)
            if xa is not None and xb is not None:
                return (xa + xb) / 2
            return xa if xa is not None else xb
        return self._ref2x.get(v.ref_pos)

    # ------------------------------------------------------------- labels

    def _update_base_labels(self) -> None:
        if self.read is None or self.read.trace is None or not len(self._xs):
            for t in self._labels:
                t.hide()
            return
        (x0, x1), _ = self.plot.getViewBox().viewRange()
        i0, i1 = np.searchsorted(self._xs, [x0, x1])
        count = i1 - i0
        if count > MAX_BASE_LABELS:
            for t in self._labels:
                t.hide()
            return
        while len(self._labels) < count:
            t = pg.TextItem(anchor=(0.5, 0.5))
            t.setFont(self._label_font)
            t.setZValue(30)
            self.plot.addItem(t)
            self._labels.append(t)
        y = self._ymax * 0.80
        k = 0
        for i in range(i0, i1):
            t = self._labels[k]
            base = self._calls[i]
            t.setText(base)
            t.setColor(base_color(base))
            t.setPos(float(self._xs[i]), y)
            t.show()
            k += 1
        for t in self._labels[k:]:
            t.hide()

    # -------------------------------------------------------------- cursor

    def center_on_ref(self, refpos: int, span_bases: int | None = None) -> None:
        x = self._nearest_x(refpos)
        if x is None:
            self.cursor_line.hide()
            return
        vb = self.plot.getViewBox()
        (x0, x1), _ = vb.viewRange()
        span = x1 - x0
        if span_bases is not None and len(self._xs) > 1:
            spacing = float(np.median(np.diff(self._xs)))
            span = span_bases * spacing
        vb.setXRange(x - span / 2, x + span / 2, padding=0)
        self.cursor_line.setPos(x)
        self.cursor_line.show()

    def set_cursor_ref(self, refpos: int) -> None:
        x = self._nearest_x(refpos)
        if x is None:
            self.cursor_line.hide()
            return
        self.cursor_line.setPos(x)
        self.cursor_line.show()
        vb = self.plot.getViewBox()
        (x0, x1), _ = vb.viewRange()
        if not (x0 <= x <= x1):
            vb.setXRange(x - (x1 - x0) / 2, x + (x1 - x0) / 2, padding=0)

    def _nearest_x(self, refpos: int) -> float | None:
        if refpos in self._ref2x:
            return self._ref2x[refpos]
        if not self._ref2x:
            return None
        aln = self.read.alignment if self.read else None
        if aln is None or not (aln.ref_start <= refpos < aln.ref_end):
            return None
        for d in range(1, 30):
            for cand in (refpos - d, refpos + d):
                if cand in self._ref2x:
                    return self._ref2x[cand]
        return None

    def _on_click(self, ev) -> None:
        if self.read is None or not len(self._xs):
            return
        vb = self.plot.getViewBox()
        if ev.button() != Qt.LeftButton:
            return
        pos = vb.mapSceneToView(ev.scenePos())
        i = int(np.searchsorted(self._xs, pos.x()))
        i = min(max(i, 0), len(self._xs) - 1)
        if i > 0 and abs(self._xs[i - 1] - pos.x()) < abs(self._xs[i] - pos.x()):
            i -= 1
        ref = int(self._refs[i]) if i < len(self._refs) else -1
        self.positionClicked.emit(ref, self.read)
