"""Main macOS desktop review window and project workflow."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, QStandardPaths, Qt, QThreadPool, QTimer, Signal
from PySide6.QtGui import QAction, QCloseEvent, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QDockWidget,
    QFileDialog,
    QInputDialog,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QSplitter,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from ..core.export import export_variants_csv
from ..core.model import Project, Read, Sample, Variant
from ..core.pipeline import load_project_inputs
from ..core.projectio import load_project, save_project
from ..core.reference import load_reference
from ..devtools.demogen import generate_demo
from .alignment_view import AlignmentStrip
from .qc_panel import QCPanel
from .sample_list import SampleList
from .summary_bar import SummaryBar
from .trace_view import TraceView
from .variant_table import VariantPanel
from .workers import AnalyzeJob


PROJECT_FILTER = "Sanger Seek projects (*.sanger-seek.json *.json)"
INPUT_FILTER = "Sanger data (*.ab1 *.abi *.ab *.seq *.fa *.fasta *.fna *.ffn *.gb *.gbk *.genbank *.gbff)"


class MainWindow(QMainWindow):
    screenshotSaved = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.project = Project()
        self.current_sample: Sample | None = None
        self.cursor_ref = 0
        self.pool = QThreadPool.globalInstance()
        # PySide QRunnables must stay strongly referenced until their queued
        # completion signal arrives; otherwise the Python wrapper can be
        # collected while Qt is still executing it.
        self._jobs: dict[tuple[int, str], AnalyzeJob] = {}
        self._pending = 0
        self._generation = 0
        self._dirty = False
        self._screenshot_path: str | None = None
        self.quit_after_screenshot = False
        self.settings = QSettings()
        self.setAcceptDrops(True)
        self.setWindowTitle("Sanger Seek")
        self.resize(1380, 900)

        self._build_ui()
        self._build_actions()
        self._restore_settings()
        self._update_ui()

    # --------------------------------------------------------------- setup

    def _build_ui(self) -> None:
        body = QWidget()
        layout = QVBoxLayout(body)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(5)
        self.summary = SummaryBar()
        self.alignment = AlignmentStrip()
        self.forward_trace = TraceView("Forward chromatogram")
        self.reverse_trace = TraceView("Reverse chromatogram")
        self.variants = VariantPanel()
        layout.addWidget(self.summary)
        layout.addWidget(self.alignment)
        traces = QSplitter(Qt.Vertical)
        traces.addWidget(self.forward_trace)
        traces.addWidget(self.reverse_trace)
        traces.setSizes([260, 260])
        lower = QSplitter(Qt.Vertical)
        lower.addWidget(traces)
        lower.addWidget(self.variants)
        lower.setSizes([550, 230])
        layout.addWidget(lower, 1)
        self.setCentralWidget(body)

        self.sample_list = SampleList()
        self.sample_list.setSelectionMode(QAbstractItemView.SingleSelection)
        self.samples_dock = QDockWidget("Samples", self)
        self.samples_dock.setObjectName("samplesDock")
        self.samples_dock.setWidget(self.sample_list)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.samples_dock)

        self.qc = QCPanel()
        self.qc_dock = QDockWidget("QC details", self)
        self.qc_dock.setObjectName("qcDock")
        self.qc_dock.setWidget(self.qc)
        self.addDockWidget(Qt.RightDockWidgetArea, self.qc_dock)
        self.qc_dock.hide()

        self.progress = QProgressBar()
        self.progress.setMaximumWidth(180)
        self.progress.hide()
        self.statusBar().addPermanentWidget(self.progress)

        self.sample_list.sampleSelected.connect(self._select_sample)
        self.summary.qcToggled.connect(self.qc_dock.setVisible)
        self.qc_dock.visibilityChanged.connect(self.summary.qc_btn.setChecked)
        self.alignment.cursorRequested.connect(self.set_cursor)
        self.forward_trace.positionClicked.connect(lambda pos, _read: self.set_cursor(pos))
        self.reverse_trace.positionClicked.connect(lambda pos, _read: self.set_cursor(pos))
        self.variants.variantSelected.connect(self._select_variant)
        self.variants.exportRequested.connect(self.export_visible)

    def _action(self, text: str, slot, shortcut=None) -> QAction:
        action = QAction(text, self)
        action.triggered.connect(slot)
        if shortcut:
            action.setShortcut(shortcut)
        return action

    def _build_actions(self) -> None:
        self.open_folder_action = self._action("Open Folder…", self.open_folder, QKeySequence.Open)
        self.add_files_action = self._action("Add Files…", self.add_files, "Ctrl+Shift+O")
        self.open_project_action = self._action("Open Project…", self.open_project, "Ctrl+Alt+O")
        self.save_action = self._action("Save Project", self.save, QKeySequence.Save)
        self.save_as_action = self._action("Save Project As…", self.save_as, QKeySequence.SaveAs)
        self.export_action = self._action("Export Visible Variants…", self.export_visible, "Ctrl+E")
        self.demo_action = self._action("Open Demo", self.open_demo)
        self.reassign_action = self._action("Move Read to Sample…", self.reassign_read)
        self.prev_base_action = self._action("Previous Base", lambda: self.move_cursor(-1), Qt.Key_Left)
        self.next_base_action = self._action("Next Base", lambda: self.move_cursor(1), Qt.Key_Right)
        self.prev_variant_action = self._action("Previous Variant", lambda: self.move_variant(-1), "Shift+Ctrl+Left")
        self.next_variant_action = self._action("Next Variant", lambda: self.move_variant(1), "Shift+Ctrl+Right")
        self.about_action = self._action("About Sanger Seek", self.about)

        file_menu = self.menuBar().addMenu("File")
        for action in (
            self.open_folder_action, self.add_files_action, self.open_project_action,
            self.save_action, self.save_as_action, self.export_action,
        ):
            file_menu.addAction(action)
        file_menu.addSeparator()
        file_menu.addAction(self.demo_action)
        file_menu.addSeparator()
        file_menu.addAction(self._action("Close Window", self.close, QKeySequence.Close))

        edit_menu = self.menuBar().addMenu("Edit")
        edit_menu.addAction(self.reassign_action)
        navigate = self.menuBar().addMenu("Navigate")
        for action in (self.prev_base_action, self.next_base_action, self.prev_variant_action, self.next_variant_action):
            navigate.addAction(action)
        view = self.menuBar().addMenu("View")
        view.addAction(self.samples_dock.toggleViewAction())
        view.addAction(self.qc_dock.toggleViewAction())
        help_menu = self.menuBar().addMenu("Help")
        help_menu.addAction(self.about_action)

        toolbar = QToolBar("Main", self)
        toolbar.setObjectName("mainToolbar")
        toolbar.setMovable(False)
        for action in (self.open_folder_action, self.add_files_action, self.save_action, self.export_action):
            toolbar.addAction(action)
        toolbar.addSeparator()
        toolbar.addAction(self.prev_variant_action)
        toolbar.addAction(self.next_variant_action)
        self.addToolBar(toolbar)

    # ----------------------------------------------------------- persistence

    def _restore_settings(self) -> None:
        geometry = self.settings.value("window/geometry")
        state = self.settings.value("window/state")
        if geometry:
            self.restoreGeometry(geometry)
        if state:
            self.restoreState(state)

    def closeEvent(self, event: QCloseEvent) -> None:
        if not self._confirm_discard():
            event.ignore()
            return
        self.settings.setValue("window/geometry", self.saveGeometry())
        self.settings.setValue("window/state", self.saveState())
        event.accept()

    def _start_dir(self) -> str:
        return str(self.settings.value("files/lastDirectory", str(Path.home())))

    def _remember_dir(self, path: str | Path) -> None:
        p = Path(path)
        self.settings.setValue("files/lastDirectory", str(p if p.is_dir() else p.parent))

    # --------------------------------------------------------------- loading

    def open_folder(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Open Sanger data folder", self._start_dir())
        if path:
            self.add_paths([path])

    def add_files(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(self, "Add Sanger data", self._start_dir(), INPUT_FILTER)
        if paths:
            self.add_paths(paths)

    def add_paths(self, paths: list[str]) -> None:
        if not paths:
            return
        self._remember_dir(paths[0])
        try:
            scan = load_project_inputs(paths, self.project)
            if scan.references:
                preferred = next((p for p in scan.references if p.suffix.lower() in {".gb", ".gbk", ".genbank", ".gbff"}), scan.references[0])
                self.project.reference = load_reference(preferred)
                if len(scan.references) > 1:
                    self.statusBar().showMessage(f"Loaded {preferred.name}; {len(scan.references) - 1} other reference file(s) ignored", 8000)
            if not self.project.samples and not self.project.reference:
                QMessageBox.information(self, "No Sanger files", "No supported Sanger or reference files were found.")
                return
            self._dirty = True
            self._begin_analysis()
        except Exception as exc:
            self._show_error("Could not load input", exc)

    def open_project(self) -> None:
        if not self._confirm_discard():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open project", self._start_dir(), PROJECT_FILTER)
        if path:
            self.open_project_path(path)

    def open_project_path(self, path: str) -> None:
        try:
            project, warnings = load_project(path)
            self.project = project
            self._dirty = False
            self._remember_dir(path)
            self._generation += 1
            self.current_sample = None
            self._begin_analysis()
            if warnings:
                QMessageBox.warning(self, "Project opened with warnings", "\n".join(warnings))
        except Exception as exc:
            self._show_error("Could not open project", exc)

    def open_demo(self) -> None:
        app_data = Path(QStandardPaths.writableLocation(QStandardPaths.AppLocalDataLocation))
        demo = app_data / "demo-data"
        try:
            demo.mkdir(parents=True, exist_ok=True)
            generate_demo(demo)
            if self.project.samples and not self._confirm_discard():
                return
            self.project = Project()
            self.add_paths([str(demo)])
        except Exception as exc:
            self._show_error("Could not create demo data", exc)

    # -------------------------------------------------------------- analysis

    def _begin_analysis(self) -> None:
        self._generation += 1
        generation = self._generation
        samples = self.project.samples
        self.sample_list.set_samples(samples, self.current_sample.key if self.current_sample else None)
        self._pending = len(samples)
        self.progress.setRange(0, max(self._pending, 1))
        self.progress.setValue(0)
        self.progress.setVisible(bool(samples))
        self.statusBar().showMessage(f"Analyzing {len(samples)} sample(s)…" if samples else "Reference loaded")
        if not samples:
            self._update_ui()
            self._maybe_screenshot()
            return
        for sample in samples:
            sample.analyzed = False
            job = AnalyzeJob(sample, self.project.reference, self.project.config)
            job.signals.finished.connect(lambda key, g=generation: self._analysis_done(key, None, g))
            job.signals.failed.connect(lambda key, error, g=generation: self._analysis_done(key, error, g))
            self._jobs[(generation, sample.key)] = job
            self.pool.start(job)

    def _analysis_done(self, key: str, error: str | None, generation: int) -> None:
        self._jobs.pop((generation, key), None)
        if generation != self._generation:
            return
        sample = self.project.sample_by_key(key)
        if sample:
            self.sample_list.refresh_sample(sample)
        self._pending = max(self._pending - 1, 0)
        self.progress.setValue(self.progress.maximum() - self._pending)
        if self.current_sample is sample:
            self._update_ui()
        if self._pending == 0:
            self.progress.hide()
            self.statusBar().showMessage("Analysis complete" if not error else f"Analysis finished with an error: {error}", 6000)
            if self.current_sample is None and self.project.samples:
                self.sample_list.setCurrentRow(0)
            self._update_ui()
            self._maybe_screenshot()

    # ------------------------------------------------------------- selection

    def _select_sample(self, key: str) -> None:
        self.current_sample = self.project.sample_by_key(key)
        if self.current_sample and self.project.reference:
            aligned = [r.alignment for r in self.current_sample.reads if r.alignment]
            if aligned and not any(a.ref_start <= self.cursor_ref < a.ref_end for a in aligned):
                self.cursor_ref = min(a.ref_start for a in aligned)
        self._update_ui()

    def _display_reads(self) -> tuple[Read | None, Read | None]:
        if self.current_sample is None:
            return None, None
        fwd = self.current_sample.forward_reads
        rev = self.current_sample.reverse_reads
        first = fwd[0] if fwd else (self.current_sample.reads[0] if self.current_sample.reads else None)
        second = rev[0] if rev else next((r for r in self.current_sample.reads if r is not first), None)
        return first, second

    def _update_ui(self) -> None:
        sample = self.current_sample
        ref = self.project.reference
        self.summary.update_summary(sample, ref)
        self.alignment.set_data(ref, sample)
        self.alignment.set_cursor(self.cursor_ref)
        first, second = self._display_reads()
        self.forward_trace.set_read(first, ref)
        self.reverse_trace.set_read(second, ref)
        variants = sample.variants if sample else []
        self.forward_trace.set_variants(variants)
        self.reverse_trace.set_variants(variants)
        self.variants.set_variants(variants)
        self.qc.set_data(sample, self.cursor_ref)
        self.set_cursor(self.cursor_ref, center=False)
        self.save_action.setEnabled(bool(self.project.samples or ref))
        self.save_as_action.setEnabled(bool(self.project.samples or ref))
        self.export_action.setEnabled(bool(variants))
        self.reassign_action.setEnabled(bool(sample and sample.reads))
        title = "Sanger Seek"
        if self.project.path:
            title += f" — {Path(self.project.path).name}"
        if self._dirty:
            title += " *"
        self.setWindowTitle(title)

    def set_cursor(self, refpos: int, center: bool = False) -> None:
        ref = self.project.reference
        if ref is None or refpos < 0:
            return
        self.cursor_ref = max(0, min(ref.n - 1, refpos))
        self.alignment.set_cursor(self.cursor_ref)
        if center:
            self.forward_trace.center_on_ref(self.cursor_ref, 45)
            self.reverse_trace.center_on_ref(self.cursor_ref, 45)
        else:
            self.forward_trace.set_cursor_ref(self.cursor_ref)
            self.reverse_trace.set_cursor_ref(self.cursor_ref)
        self.qc.set_data(self.current_sample, self.cursor_ref)

    def move_cursor(self, delta: int) -> None:
        self.set_cursor(self.cursor_ref + delta)

    def _select_variant(self, variant: Variant) -> None:
        self.set_cursor(variant.ref_pos, center=True)

    def move_variant(self, delta: int) -> None:
        variants = self.variants.visible_variants()
        if not variants:
            return
        positions = [v.ref_pos for v in variants]
        if delta > 0:
            chosen = next((v for v in variants if v.ref_pos > self.cursor_ref), variants[0])
        else:
            chosen = next((v for v in reversed(variants) if v.ref_pos < self.cursor_ref), variants[-1])
        self.variants.select_variant(chosen)
        self._select_variant(chosen)

    # --------------------------------------------------------------- project

    def save(self) -> None:
        if not self.project.path:
            self.save_as()
            return
        try:
            save_project(self.project, self.project.path)
            self._dirty = False
            self._update_ui()
            self.statusBar().showMessage("Project saved", 4000)
        except Exception as exc:
            self._show_error("Could not save project", exc)

    def save_as(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save project", self._start_dir(), PROJECT_FILTER)
        if not path:
            return
        if not path.endswith(".json"):
            path += ".sanger-seek.json"
        self.project.path = path
        self._remember_dir(path)
        self.save()

    def export_visible(self) -> None:
        if self.current_sample is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export visible variants", self._start_dir(), "CSV files (*.csv)")
        if not path:
            return
        if not path.lower().endswith(".csv"):
            path += ".csv"
        visible_ids = {v.id for v in self.variants.visible_variants()}
        try:
            count = export_variants_csv(
                path, [self.current_sample], self.project.reference,
                variants_filter=lambda v: v.id in visible_ids,
            )
            self._remember_dir(path)
            self.statusBar().showMessage(f"Exported {count} variant(s) to {Path(path).name}", 6000)
        except Exception as exc:
            self._show_error("Could not export variants", exc)

    def reassign_read(self) -> None:
        source = self.current_sample
        if source is None or not source.reads:
            return
        labels = [r.label for r in source.reads]
        label, ok = QInputDialog.getItem(self, "Move read", "Read:", labels, 0, False)
        if not ok:
            return
        targets = [s for s in self.project.samples if s is not source]
        names = [s.name for s in targets] + ["New sample…"]
        target_name, ok = QInputDialog.getItem(self, "Move read", "Destination sample:", names, 0, False)
        if not ok:
            return
        read = next(r for r in source.reads if r.label == label)
        if target_name == "New sample…":
            name, ok = QInputDialog.getText(self, "New sample", "Sample name:")
            if not ok or not name.strip():
                return
            key = name.strip()
            if self.project.sample_by_key(key):
                QMessageBox.warning(self, "Sample exists", "A sample with that key already exists.")
                return
            target = Sample(key=key, name=key)
            self.project.samples.append(target)
        else:
            target = targets[names.index(target_name)]
        source.reads.remove(read)
        read.id = f"{target.key}/{read.label.lower()}"
        target.reads.append(read)
        if not source.reads:
            self.project.samples.remove(source)
        self.current_sample = target
        self._dirty = True
        self.project.samples.sort(key=lambda s: s.name.lower())
        self._begin_analysis()

    # ----------------------------------------------------------- mac / misc

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls() and all(url.isLocalFile() for url in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:
        paths = [url.toLocalFile() for url in event.mimeData().urls() if url.isLocalFile()]
        if len(paths) == 1 and paths[0].endswith(".json"):
            if self._confirm_discard():
                self.open_project_path(paths[0])
        else:
            self.add_paths(paths)
        event.acceptProposedAction()

    def open_external_path(self, path: str) -> None:
        """Open a path delivered by Finder or another macOS application."""
        if path.endswith(".json"):
            if self._confirm_discard():
                self.open_project_path(path)
        else:
            self.add_paths([path])

    def request_screenshot(self, path: str) -> None:
        self._screenshot_path = str(Path(path).expanduser().resolve())
        if self._pending == 0:
            QTimer.singleShot(400, self._maybe_screenshot)

    def _maybe_screenshot(self) -> None:
        if not self._screenshot_path or self._pending:
            return
        path = Path(self._screenshot_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.grab().save(str(path), "PNG")
        self.statusBar().showMessage(f"Screenshot saved to {path}", 6000)
        self._screenshot_path = None
        self.screenshotSaved.emit(str(path))
        if self.quit_after_screenshot:
            # Completion signals are queued just before QRunnable.run returns;
            # let Qt finish releasing the pooled workers before tearing down
            # the application in automated screenshot mode.
            QTimer.singleShot(750, self._finish_screenshot_mode)

    def _finish_screenshot_mode(self) -> None:
        self._dirty = False
        self.close()
        QApplication.instance().quit()

    def _confirm_discard(self) -> bool:
        if not self._dirty:
            return True
        answer = QMessageBox.question(
            self, "Unsaved project", "Discard unsaved project changes?",
            QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Cancel,
        )
        return answer == QMessageBox.Discard

    def _show_error(self, title: str, exc: Exception) -> None:
        QMessageBox.critical(self, title, str(exc))
        self.statusBar().showMessage(str(exc), 8000)

    def about(self) -> None:
        QMessageBox.about(
            self,
            "About Sanger Seek",
            "<b>Sanger Seek 0.1.0</b><br><br>Local Sanger chromatogram review, "
            "strand reconciliation, and candidate variant annotation.<br><br>"
            "Sequence data remains on this Mac. Results are for research review and are not "
            "a clinically validated diagnostic report.",
        )
