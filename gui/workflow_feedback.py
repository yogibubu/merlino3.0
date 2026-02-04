"""gui/workflow_feedback.py

Modal report viewer for interactive workflow confirmation.

Behavior:
- ask_open_report(): asks user if they want to inspect a report file
- if YES: opens a modal dialog with a scrollable QTextEdit and OK / Not OK buttons
- returns True if OK, False otherwise (including if user chooses not to open the report)
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtGui import QFont, QTextCursor
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
)


class ReportViewer(QDialog):
    """Modal, scrollable, read-only report viewer with explicit OK / Not OK."""

    def __init__(self, title: str, path: Path, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(800, 650)

        layout = QVBoxLayout(self)

        self.text = QTextEdit(self)
        self.text.setReadOnly(True)
        mono = QFont("Courier New")
        mono.setStyleHint(QFont.Monospace)
        self.text.setFont(mono)
        layout.addWidget(self.text)

        if path.exists():
            text = path.read_text(encoding="utf-8", errors="replace")
        else:
            text = f"Report not found: {path}"

        self.text.setPlainText(text)
        self.text.moveCursor(QTextCursor.Start)
        self.text.setFocus()

        btns = QHBoxLayout()
        self.btn_ok = QPushButton("OK", self)
        self.btn_no = QPushButton("Not OK", self)
        btns.addWidget(self.btn_ok)
        btns.addWidget(self.btn_no)
        layout.addLayout(btns)

        self.btn_ok.clicked.connect(self.accept)
        self.btn_no.clicked.connect(self.reject)


def ask_open_report(
    parent,
    title: str,
    ok_message: str,
    report_filename: str,
    working_dir: Path | None = None,
) -> bool:
    """Ask the user if they want to inspect a report.

    Returns
    -------
    bool
        True if user opened the report and clicked OK.
        False otherwise.
    """

    reply = QMessageBox.question(
        parent,
        title,
        f"{ok_message}\n\nDo you want to inspect the report?",
        QMessageBox.Yes | QMessageBox.No,
        QMessageBox.Yes,
    )

    if reply != QMessageBox.Yes:
        return True  # user trusts it; continue workflow

    base = working_dir if working_dir is not None else (Path.cwd() / "working")
    try:
        base.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    path = base / report_filename

    viewer = ReportViewer(title, path, parent=parent)
    res = viewer.exec()
    return (res == QDialog.Accepted)
