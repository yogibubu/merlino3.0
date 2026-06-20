from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QFileDialog,
    QCheckBox,
    QComboBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from advanced.launchers.puckering_dvr_launcher import PuckeringDVRLauncher


class DVRWindow(QMainWindow):
    """Dedicated GUI for DVR analysis of completed Gaussian scan/path logs."""

    def __init__(self, workdir: Path, repo_root: Path | None = None, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.repo_root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]
        self.launcher: PuckeringDVRLauncher | None = None

        self.setWindowTitle("Merlino DVR - Gaussian Scan Analysis")
        self.resize(760, 520)
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        header = QLabel("Path DVR - Gaussian Scan Analysis")
        header.setStyleSheet("font-size: 18px; font-weight: bold;")
        layout.addWidget(header)

        note = QLabel(
            "The DVR backend reads completed Gaussian scan/path logs. "
            "Gaussian input/path generation is handled by Merlino, not by DVR."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        input_group = QGroupBox("Input")
        input_layout = QVBoxLayout(input_group)

        row_log = QHBoxLayout()
        row_log.addWidget(QLabel("Gaussian log:"))
        self.log_edit = QLineEdit(str(self.workdir / "gauin.log"))
        row_log.addWidget(self.log_edit)
        browse_log = QPushButton("Browse")
        browse_log.clicked.connect(self._browse_log)
        row_log.addWidget(browse_log)
        input_layout.addLayout(row_log)

        layout.addWidget(input_group)

        output_group = QGroupBox("Outputs")
        output_layout = QVBoxLayout(output_group)

        row_out = QHBoxLayout()
        row_out.addWidget(QLabel("Output dir:"))
        self.outdir_edit = QLineEdit(str(self.workdir / "puckering_dvr_outputs"))
        row_out.addWidget(self.outdir_edit)
        browse_out = QPushButton("Browse")
        browse_out.clicked.connect(self._browse_outdir)
        row_out.addWidget(browse_out)
        output_layout.addLayout(row_out)

        row_fig = QHBoxLayout()
        row_fig.addWidget(QLabel("Figure dir:"))
        self.figdir_edit = QLineEdit(str(self.workdir / "puckering_dvr_figs"))
        row_fig.addWidget(self.figdir_edit)
        browse_fig = QPushButton("Browse")
        browse_fig.clicked.connect(self._browse_figdir)
        row_fig.addWidget(browse_fig)
        output_layout.addLayout(row_fig)

        row_prefix = QHBoxLayout()
        row_prefix.addWidget(QLabel("Prefix:"))
        self.prefix_edit = QLineEdit("puckering_dvr")
        row_prefix.addWidget(self.prefix_edit)
        output_layout.addLayout(row_prefix)

        layout.addWidget(output_group)

        settings_group = QGroupBox("Hamiltonian Settings")
        settings_layout = QVBoxLayout(settings_group)

        row_method = QHBoxLayout()
        row_method.addWidget(QLabel("Boundary:"))
        self.boundary_combo = QComboBox()
        self.boundary_combo.addItems(["periodic", "nonperiodic"])
        row_method.addWidget(self.boundary_combo)

        row_method.addWidget(QLabel("Solver:"))
        self.solver_combo = QComboBox()
        self.solver_combo.addItems(["auto", "fourier", "gaussian", "sinc-dvr"])
        self.solver_combo.setCurrentText("fourier")
        row_method.addWidget(self.solver_combo)
        settings_layout.addLayout(row_method)

        row_flags = QHBoxLayout()
        self.rotconst_check = QCheckBox("Compute rotational constants")
        self.rotconst_check.setChecked(True)
        row_flags.addWidget(self.rotconst_check)

        self.cremer_check = QCheckBox("Label Cremer-Pople")
        self.cremer_check.setChecked(False)
        row_flags.addWidget(self.cremer_check)
        settings_layout.addLayout(row_flags)

        layout.addWidget(settings_group)

        actions = QHBoxLayout()
        self.run_button = QPushButton("Run DVR")
        self.run_button.clicked.connect(self.run_dvr)
        actions.addWidget(self.run_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        actions.addWidget(close_button)
        layout.addLayout(actions)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setPlaceholderText("DVR command output will appear here.")
        layout.addWidget(self.output_text)

    def run_dvr(self) -> None:
        log_path = Path(self.log_edit.text().strip()).expanduser()
        outdir = Path(self.outdir_edit.text().strip()).expanduser()
        figdir = Path(self.figdir_edit.text().strip()).expanduser()
        prefix = self.prefix_edit.text().strip() or "puckering_dvr"

        outdir.mkdir(parents=True, exist_ok=True)
        figdir.mkdir(parents=True, exist_ok=True)

        self.run_button.setEnabled(False)
        self.output_text.setPlainText("Starting DVR analysis...")

        self.launcher = PuckeringDVRLauncher(self.workdir, self.repo_root, parent=self)
        self.launcher.finished.connect(self._on_finished)
        self.launcher.start_path_analysis(
            log_path,
            outdir,
            figdir,
            prefix,
            self.boundary_combo.currentText(),
            self.solver_combo.currentText(),
            compute_rotconst=self.rotconst_check.isChecked(),
            label_cremer_pople=self.cremer_check.isChecked(),
        )

    def _on_finished(self, success: bool, message: str) -> None:
        self.run_button.setEnabled(True)
        self.output_text.setPlainText(message)
        if success:
            QMessageBox.information(self, "Path DVR", message)
        else:
            QMessageBox.critical(self, "Path DVR failed", message)

    def _browse_log(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Gaussian log", str(self.workdir), "Gaussian Logs (*.log *.out);;All Files (*)"
        )
        if path:
            self.log_edit.setText(path)

    def _browse_outdir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select DVR output directory", str(self.workdir))
        if path:
            self.outdir_edit.setText(path)

    def _browse_figdir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select DVR figure directory", str(self.workdir))
        if path:
            self.figdir_edit.setText(path)

