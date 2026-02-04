from pathlib import Path

from PySide6.QtWidgets import (
    QMainWindow, QTextEdit, QTabWidget,
    QWidget, QVBoxLayout, QPushButton
)


class ProvaViewer(QMainWindow):
    """
    Full-window viewer for prova outputs:
    - provout
    - gauin.gjf

    Explicitly launches Gaussian via callback.
    """

    def __init__(self, workdir: Path, run_gaussian_callback, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self._run_gaussian_callback = run_gaussian_callback

        self.setWindowTitle("Prova Viewer")
        self.resize(1100, 850)

        # ----------------------------------------------------------
        # Central widget & layout
        # ----------------------------------------------------------
        central = QWidget()
        layout = QVBoxLayout(central)
        self.setCentralWidget(central)

        # ----------------------------------------------------------
        # Tabs
        # ----------------------------------------------------------
        tabs = QTabWidget()
        layout.addWidget(tabs)

        # provout
        provout_edit = QTextEdit()
        provout_edit.setReadOnly(True)
        provout = self.workdir / "provout"
        if provout.exists():
            provout_edit.setPlainText(provout.read_text())
        tabs.addTab(provout_edit, "provout")

        # gauin.gjf
        gauin_edit = QTextEdit()
        gauin_edit.setReadOnly(True)
        gauin = self.workdir / "gauin.gjf"
        if gauin.exists():
            gauin_edit.setPlainText(gauin.read_text())
        tabs.addTab(gauin_edit, "gauin.gjf")

        # ----------------------------------------------------------
        # Run Gaussian button
        # ----------------------------------------------------------
        run_btn = QPushButton("Run Gaussian")
        run_btn.setToolTip("Launch Gaussian (gdv) using gauin.gjf")
        run_btn.clicked.connect(self._on_run_gaussian)
        layout.addWidget(run_btn)

    # --------------------------------------------------------------
    def _on_run_gaussian(self):
        """
        Explicit user action: launch Gaussian and close viewer.
        """
        self._run_gaussian_callback()
        self.close()

