from __future__ import annotations

from datetime import datetime
import hashlib
import json
from pathlib import Path
import sys

from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices
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

from advanced.launchers.gaussian_launcher import GaussianLauncher
from advanced.launchers.puckering_dvr_launcher import PuckeringDVRLauncher


class DVRWindow(QMainWindow):
    """Dedicated GUI for DVR analysis of completed Gaussian scan/path logs."""

    def __init__(self, workdir: Path, repo_root: Path | None = None, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.repo_root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]
        self.launcher: PuckeringDVRLauncher | None = None
        self.gaussian_launcher: GaussianLauncher | None = None
        self._run_dvr_after_gaussian = False

        self.setWindowTitle("Merlino DVR - Gaussian Scan Analysis")
        self.resize(900, 760)
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

        row_log_actions = QHBoxLayout()
        latest_log = QPushButton("Use Latest Gaussian Log")
        latest_log.clicked.connect(lambda: self._select_latest_log())
        row_log_actions.addWidget(latest_log)
        preview_log = QPushButton("Preflight / Preview")
        preview_log.clicked.connect(lambda: self.preview_log())
        row_log_actions.addWidget(preview_log)
        row_log_actions.addStretch()
        input_layout.addLayout(row_log_actions)

        self.preflight_label = QLabel("Preflight not run.")
        self.preflight_label.setWordWrap(True)
        input_layout.addWidget(self.preflight_label)

        layout.addWidget(input_group)

        workflow_group = QGroupBox("Guided Workflow")
        workflow_layout = QVBoxLayout(workflow_group)
        workflow_note = QLabel(
            "Use Merlino to prepare gauin.gjf, then run Gaussian here and pass the completed log to DVR. "
            "DVR only reads Gaussian output."
        )
        workflow_note.setWordWrap(True)
        workflow_layout.addWidget(workflow_note)

        row_workflow = QHBoxLayout()
        self.gaussian_button = QPushButton("Run Gaussian")
        self.gaussian_button.clicked.connect(self.run_gaussian)
        row_workflow.addWidget(self.gaussian_button)
        self.gaussian_dvr_button = QPushButton("Run Gaussian Then DVR")
        self.gaussian_dvr_button.clicked.connect(self.run_gaussian_then_dvr)
        row_workflow.addWidget(self.gaussian_dvr_button)
        row_workflow.addStretch()
        workflow_layout.addLayout(row_workflow)

        layout.addWidget(workflow_group)

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
        self.solver_combo.addItems([
            "auto",
            "fourier",
            "gaussian",
            "sinc-dvr",
            "fortran-sinc-dvr",
            "fortran-gaussian",
        ])
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

        results_group = QGroupBox("Results")
        results_layout = QHBoxLayout(results_group)
        refresh_results = QPushButton("Refresh Results")
        refresh_results.clicked.connect(self.refresh_results)
        results_layout.addWidget(refresh_results)
        open_summary = QPushButton("Open Summary")
        open_summary.clicked.connect(lambda: self._open_result("summary"))
        results_layout.addWidget(open_summary)
        open_levels = QPushButton("Open Levels")
        open_levels.clicked.connect(lambda: self._open_result("levels"))
        results_layout.addWidget(open_levels)
        open_profile = QPushButton("Open Profile")
        open_profile.clicked.connect(lambda: self._open_result("profile"))
        results_layout.addWidget(open_profile)
        open_outdir = QPushButton("Open Output Dir")
        open_outdir.clicked.connect(lambda: self._open_path(Path(self.outdir_edit.text().strip()).expanduser()))
        results_layout.addWidget(open_outdir)
        open_figdir = QPushButton("Open Figure Dir")
        open_figdir.clicked.connect(lambda: self._open_path(Path(self.figdir_edit.text().strip()).expanduser()))
        results_layout.addWidget(open_figdir)
        layout.addWidget(results_group)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setPlaceholderText("DVR command output will appear here.")
        layout.addWidget(self.output_text)

    def run_dvr(self) -> None:
        log_path = Path(self.log_edit.text().strip()).expanduser()
        outdir = Path(self.outdir_edit.text().strip()).expanduser()
        figdir = Path(self.figdir_edit.text().strip()).expanduser()
        prefix = self.prefix_edit.text().strip() or "puckering_dvr"
        if not log_path.exists():
            QMessageBox.critical(self, "Path DVR failed", f"Gaussian log not found: {log_path}")
            return

        outdir.mkdir(parents=True, exist_ok=True)
        figdir.mkdir(parents=True, exist_ok=True)

        self.run_button.setEnabled(False)

        self.launcher = PuckeringDVRLauncher(self.workdir, self.repo_root, parent=self)
        self.launcher.finished.connect(self._on_finished)
        args = self.launcher.build_path_analysis_args(
            log_path,
            outdir,
            figdir,
            prefix,
            self.boundary_combo.currentText(),
            self.solver_combo.currentText(),
            compute_rotconst=self.rotconst_check.isChecked(),
            label_cremer_pople=self.cremer_check.isChecked(),
        )
        manifest = self._write_run_manifest(log_path, outdir, figdir, prefix, args)
        self.output_text.setPlainText(
            "Starting DVR analysis...\n\n"
            f"Command: {sys.executable} {' '.join(args)}\n"
            f"Manifest: {manifest}\n"
        )
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
        if success:
            self.refresh_results(prefix_message=message)
        else:
            self.output_text.setPlainText(message)
        self._run_dvr_after_gaussian = False
        if success:
            QMessageBox.information(self, "Path DVR", message)
        else:
            QMessageBox.critical(self, "Path DVR failed", message)

    def run_gaussian(self) -> None:
        self._run_dvr_after_gaussian = False
        self._start_gaussian()

    def run_gaussian_then_dvr(self) -> None:
        self._run_dvr_after_gaussian = True
        self._start_gaussian()

    def _start_gaussian(self) -> None:
        self.gaussian_button.setEnabled(False)
        self.gaussian_dvr_button.setEnabled(False)
        self.output_text.setPlainText("Starting Gaussian from Merlino workdir...\n")
        self.gaussian_launcher = GaussianLauncher(self.workdir, parent=self)
        self.gaussian_launcher.finished.connect(self._on_gaussian_finished)
        self.gaussian_launcher.start()

    def _on_gaussian_finished(self, success: bool, message: str) -> None:
        self.gaussian_button.setEnabled(True)
        self.gaussian_dvr_button.setEnabled(True)
        self.output_text.setPlainText(message)
        if success:
            self._select_latest_log(show_message=False)
            self.preview_log(show_message=False)
            if self._run_dvr_after_gaussian:
                self.run_dvr()
        else:
            self._run_dvr_after_gaussian = False
            QMessageBox.critical(self, "Gaussian failed", message)

    def preview_log(self, show_message: bool = True) -> None:
        log_path = Path(self.log_edit.text().strip()).expanduser()
        ok, lines = self._scan_log_summary(log_path)
        text = "\n".join(lines)
        self.preflight_label.setText(lines[0] if lines else "Preflight failed.")
        self.output_text.setPlainText(text)
        if show_message and not ok:
            QMessageBox.warning(self, "Gaussian log preflight", text)

    def refresh_results(self, prefix_message: str | None = None) -> None:
        paths = self._result_paths()
        chunks = []
        if prefix_message:
            chunks.append(prefix_message.strip())
        chunks.append("DVR Results")
        for label in ("summary", "levels", "profile", "manifest"):
            path = paths.get(label)
            if path:
                chunks.append(f"{label}: {path}")
        figures = paths.get("figures", [])
        if figures:
            chunks.append("figures:")
            chunks.extend(str(path) for path in figures[:8])

        summary = paths.get("summary")
        if summary:
            chunks.append("\nSummary preview:")
            chunks.append(self._read_preview(summary, max_lines=40))
        levels = paths.get("levels")
        if levels:
            chunks.append("\nLevels preview:")
            chunks.append(self._read_preview(levels, max_lines=16))
        if not summary and not levels:
            chunks.append("No DVR result files found for the current prefix.")
        self.output_text.setPlainText("\n".join(chunks))

    def _select_latest_log(self, show_message: bool = True) -> None:
        candidates = self._candidate_logs()
        if not candidates:
            msg = f"No Gaussian .log/.out files found in {self.workdir}"
            self.preflight_label.setText(msg)
            if show_message:
                QMessageBox.warning(self, "Gaussian log", msg)
            return
        latest = max(candidates, key=lambda p: (p.stat().st_mtime, p.stat().st_size))
        self.log_edit.setText(str(latest))
        self.preflight_label.setText(f"Selected latest log: {latest.name}")

    def _candidate_logs(self) -> list[Path]:
        seen: set[Path] = set()
        candidates: list[Path] = []
        for name in ("gauin.log", "gauout.log"):
            path = self.workdir / name
            if path.exists():
                seen.add(path.resolve())
                candidates.append(path)
        for pattern in ("*.log", "*.out"):
            for path in self.workdir.glob(pattern):
                resolved = path.resolve()
                if path.is_file() and resolved not in seen:
                    seen.add(resolved)
                    candidates.append(path)
        return candidates

    def _scan_log_summary(self, log_path: Path) -> tuple[bool, list[str]]:
        if not log_path.exists():
            return False, [f"Missing Gaussian log: {log_path}"]
        text = log_path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        orientation_count = sum(
            1 for line in lines if "Input orientation:" in line or "Standard orientation:" in line
        )
        scf_count = sum(1 for line in lines if "SCF Done:" in line)
        normal = "Normal termination of Gaussian" in text
        scan_markers = sum(
            1 for line in lines if "Step number" in line or "Optimization completed" in line or "Scan" in line
        )
        gic_markers = sum(1 for line in lines if "QPck" in line or "PhiP" in line or "RPck" in line)
        ok = orientation_count > 0 and (scf_count > 0 or scan_markers > 0)
        status = "OK" if ok else "WARNING"
        return ok, [
            f"{status}: {log_path}",
            f"size: {log_path.stat().st_size} bytes",
            f"geometries/orientations: {orientation_count}",
            f"SCF energies: {scf_count}",
            f"scan markers: {scan_markers}",
            f"puckering GIC markers: {gic_markers}",
            f"normal termination: {'yes' if normal else 'no'}",
        ]

    def _write_run_manifest(
        self,
        log_path: Path,
        outdir: Path,
        figdir: Path,
        prefix: str,
        args: list[str],
    ) -> Path:
        manifest = {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "workdir": str(self.workdir),
            "repo_root": str(self.repo_root),
            "python_executable": sys.executable,
            "command": [sys.executable, *args],
            "gaussian_log": str(log_path),
            "gaussian_log_sha256": self._sha256(log_path),
            "outdir": str(outdir),
            "figdir": str(figdir),
            "prefix": prefix,
            "boundary": self.boundary_combo.currentText(),
            "solver": self.solver_combo.currentText(),
            "compute_rotconst": self.rotconst_check.isChecked(),
            "label_cremer_pople": self.cremer_check.isChecked(),
        }
        path = outdir / f"{prefix}_dvr_run_manifest.json"
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return path

    def _result_paths(self) -> dict[str, Path | list[Path]]:
        outdir = Path(self.outdir_edit.text().strip()).expanduser()
        figdir = Path(self.figdir_edit.text().strip()).expanduser()
        prefix = self.prefix_edit.text().strip() or "puckering_dvr"

        def first(patterns: list[str]) -> Path | None:
            matches: list[Path] = []
            for pattern in patterns:
                matches.extend(outdir.glob(pattern))
            return max(matches, key=lambda p: p.stat().st_mtime) if matches else None

        result: dict[str, Path | list[Path]] = {}
        summary = first([f"{prefix}_summary.txt", f"{prefix}_*_summary.txt"])
        levels = first([f"{prefix}_levels.csv", f"{prefix}_*_levels.csv"])
        profile = first([f"{prefix}_model_profile.csv", f"{prefix}_profile.csv", f"{prefix}_*_grid.csv"])
        manifest = first([f"{prefix}_dvr_run_manifest.json"])
        figures = sorted(
            [*figdir.glob(f"{prefix}*.pdf"), *figdir.glob(f"{prefix}*.png")],
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if summary:
            result["summary"] = summary
        if levels:
            result["levels"] = levels
        if profile:
            result["profile"] = profile
        if manifest:
            result["manifest"] = manifest
        if figures:
            result["figures"] = figures
        return result

    def _open_result(self, label: str) -> None:
        path = self._result_paths().get(label)
        if isinstance(path, Path):
            self._open_path(path)
        else:
            QMessageBox.warning(self, "DVR results", f"No {label} file found for current prefix.")

    def _open_path(self, path: Path) -> None:
        if not path.exists():
            QMessageBox.warning(self, "Open path", f"Path not found: {path}")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    @staticmethod
    def _read_preview(path: Path, max_lines: int) -> str:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        shown = lines[:max_lines]
        if len(lines) > max_lines:
            shown.append(f"... ({len(lines) - max_lines} more lines)")
        return "\n".join(shown)

    @staticmethod
    def _sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

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
