from pathlib import Path
import subprocess

from PySide6.QtWidgets import (
    QMainWindow, QTextEdit, QTabWidget, QMessageBox
)


class GaussianViewer(QMainWindow):
    """
    Full-window viewer for Gaussian outputs:
    - gauout.log
    - gicforge.fchk/prova.fchk (generated via formchk)
    """

    def __init__(self, workdir: Path, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)

        self.setWindowTitle("Gaussian Viewer")
        self.resize(1100, 850)

        # ---------------- formchk ----------------
        chk = self.workdir / "gicforge.chk"
        fchk = self.workdir / "gicforge.fchk"
        if not chk.exists():
            chk = self.workdir / "prova.chk"
            fchk = self.workdir / "prova.fchk"

        if chk.exists() and not fchk.exists():
            try:
                subprocess.run(
                    ["formchk", str(chk), str(fchk)],
                    cwd=self.workdir,
                    check=True,
                )
            except Exception as e:
                QMessageBox.warning(
                    self,
                    "formchk failed",
                    f"Could not generate {fchk.name}:\n{e}",
                )

        tabs = QTabWidget()
        self.setCentralWidget(tabs)

        # ---------------- Gaussian logs ----------------
        for log_name in ("gauin.log", "gauout.log"):
            edit = QTextEdit()
            edit.setReadOnly(True)
            log_path = self.workdir / log_name
            if log_path.exists():
                edit.setPlainText(log_path.read_text(errors="replace"))
            else:
                edit.setPlainText(f"{log_name} not found")
            tabs.addTab(edit, log_name)

        # ---------------- formatted checkpoint ----------------
        fchk_edit = QTextEdit()
        fchk_edit.setReadOnly(True)
        if fchk.exists():
            fchk_edit.setPlainText(fchk.read_text())
        tabs.addTab(fchk_edit, fchk.name)
