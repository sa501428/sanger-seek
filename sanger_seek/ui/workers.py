"""Background analysis jobs (QThreadPool) so the UI never blocks."""

from __future__ import annotations

import traceback

from PySide6.QtCore import QObject, QRunnable, Signal

from ..core.model import Config, Reference, Sample
from ..core.pipeline import analyze_sample


class AnalyzeSignals(QObject):
    finished = Signal(str)          # sample key
    failed = Signal(str, str)       # sample key, error


class AnalyzeJob(QRunnable):
    """Analyzes one sample. The sample object is owned by the job until the
    finished/failed signal is delivered on the UI thread."""

    def __init__(self, sample: Sample, reference: Reference | None, config: Config):
        super().__init__()
        self.sample = sample
        self.reference = reference
        self.config = config
        self.signals = AnalyzeSignals()

    def run(self) -> None:  # executed on a worker thread
        try:
            analyze_sample(self.sample, self.reference, self.config)
            self.signals.finished.emit(self.sample.key)
        except Exception as e:
            traceback.print_exc()
            self.sample.error = str(e)
            self.signals.failed.emit(self.sample.key, str(e))
