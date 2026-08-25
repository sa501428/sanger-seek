"""QApplication bootstrap and command-line entry point."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from PySide6.QtCore import QCoreApplication, QEvent, Qt, QTimer, Signal
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication

from .ui.main_window import MainWindow
from .ui.theme import APP_QSS


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Review Sanger chromatograms locally")
    parser.add_argument("paths", nargs="*", help="AB1/SEQ/reference files or folders to open")
    parser.add_argument("--project", help="open a .sanger-seek.json project")
    parser.add_argument("--demo", action="store_true", help="generate and open the local demo dataset")
    parser.add_argument(
        "--screenshot",
        metavar="PNG",
        help="save a screenshot after loading finishes (primarily for release QA)",
    )
    return parser


def resource_path(*parts: str) -> Path:
    return Path(__file__).resolve().parent.joinpath("resources", *parts)


class SangerSeekApplication(QApplication):
    """Capture Finder open-file events used by a native macOS bundle."""

    fileOpened = Signal(str)

    def event(self, event) -> bool:
        if event.type() == QEvent.FileOpen and event.file():
            self.fileOpened.emit(event.file())
            return True
        return super().event(event)


def create_application(argv: list[str] | None = None) -> QApplication:
    QCoreApplication.setOrganizationName("Sanger Seek")
    QCoreApplication.setOrganizationDomain("sanger-seek.local")
    QCoreApplication.setApplicationName("Sanger Seek")
    QCoreApplication.setApplicationVersion("0.1.0")
    QApplication.setAttribute(Qt.AA_DontShowIconsInMenus, False)
    app = SangerSeekApplication(argv if argv is not None else sys.argv)
    app.setApplicationDisplayName("Sanger Seek")
    app.setStyleSheet(APP_QSS)
    icon = resource_path("app-icon.png")
    if icon.exists():
        app.setWindowIcon(QIcon(str(icon)))
    return app


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    app = create_application([sys.argv[0], *(argv or [])] if argv is not None else None)
    window = MainWindow()
    app.fileOpened.connect(lambda path: window.open_external_path(path))
    if args.screenshot:
        window.quit_after_screenshot = True
    window.show()

    def load_initial() -> None:
        if args.project:
            window.open_project_path(args.project)
        elif args.demo:
            window.open_demo()
        elif args.paths:
            window.add_paths(args.paths)
        if args.screenshot:
            window.request_screenshot(args.screenshot)

    QTimer.singleShot(0, load_initial)
    return app.exec()
