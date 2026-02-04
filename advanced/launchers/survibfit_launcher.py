from __future__ import annotations

from pathlib import Path
import os
from PySide6.QtCore import QObject, QProcess, Signal


class SurvibfitLauncher(QObject):
    """
    Asynchronous Survibfit launcher using QProcess.
    """

    finished = Signal(bool, str)  # success, message

    def __init__(self, workdir: Path, merlino_fit_root: Path, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.merlino_fit_root = Path(merlino_fit_root)
        self.process = QProcess(self)

    def start_vib(self, log_path: Path, fchk_path: Path | None, out_prefix: str, scale_json: Path | None):
        if not log_path.exists():
            self.finished.emit(False, f"Log not found: {log_path}")
            return

        executable = "python"
        args = [
            "-m",
            "survibfit.cli",
            "vib",
            "--log",
            str(log_path),
            "--out",
            str(self.workdir / out_prefix),
        ]
        if fchk_path:
            args += ["--fchk", str(fchk_path)]
        if scale_json:
            args += ["--scale-json", str(scale_json)]

        self._start_process(executable, args)

    def start_gic(self, xyz_path: Path, out_path: Path, include_frag: bool = False):
        if not xyz_path.exists():
            self.finished.emit(False, f"XYZ not found: {xyz_path}")
            return

        executable = "python"
        args = [
            "-m",
            "survibfit.cli",
            "gic",
            "--xyz",
            str(xyz_path),
            "--out",
            str(out_path),
        ]
        if include_frag:
            args.append("--include-frag")

        self._start_process(executable, args)

    def _start_process(self, executable, args):
        self.process.setWorkingDirectory(str(self.workdir))
        env = QProcess.systemEnvironment()
        # prepend merlino_fit to PYTHONPATH
        py_path = str(self.merlino_fit_root)
        env = [e for e in env if not e.startswith("PYTHONPATH=")]
        env.append("PYTHONPATH=" + py_path)
        self.process.setEnvironment(env)

        self.process.finished.connect(self._on_finished)
        self.process.errorOccurred.connect(self._on_error)
        self.process.start(executable, args)

        if not self.process.waitForStarted(3000):
            self.finished.emit(False, "Failed to start Survibfit process")

    def _on_finished(self, exit_code, exit_status):
        if exit_code == 0:
            self.finished.emit(True, "Survibfit vibrational analysis completed")
        else:
            self.finished.emit(False, f"Survibfit failed (exit_code={exit_code})")

    def _on_error(self, error):
        self.finished.emit(False, f"Survibfit process error: {error}")
