from pathlib import Path

from PySide6.QtWidgets import (
    QMainWindow, QTextEdit, QTabWidget,
    QWidget, QVBoxLayout, QPushButton
)


class GICForgeViewer(QMainWindow):
    """
    Full-window viewer for GICForge outputs:
    - gicforge.out/provout
    - gauin.gjf

    Explicitly launches Gaussian via callback.
    """

    def __init__(self, workdir: Path, run_gaussian_callback, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self._run_gaussian_callback = run_gaussian_callback

        self.setWindowTitle("GICForge Viewer")
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

        # readable GICForge output
        output_edit = QTextEdit()
        output_edit.setReadOnly(True)
        output = self.workdir / "gicforge.out"
        if not output.exists():
            output = self.workdir / "provout"
        if output.exists():
            output_edit.setPlainText(output.read_text())
        tabs.addTab(output_edit, output.name if output.exists() else "gicforge.out")

        # optional B matrix / geometry-update output
        bmat_edit = QTextEdit()
        bmat_edit.setReadOnly(True)
        bmat = self.workdir / "bmat.out"
        if bmat.exists():
            bmat_edit.setPlainText(bmat.read_text())
        tabs.addTab(bmat_edit, "bmat.out")

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
