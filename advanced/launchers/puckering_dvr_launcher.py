from __future__ import annotations

from pathlib import Path
import sys

from PySide6.QtCore import QObject, QProcess, Signal
from merlino_dvr import (
    DVRRequest,
    build_fortran_bridge_args,
    build_fortran_shell_command,
    build_path_analysis_args,
    is_fortran_solver,
    resolve_dvr_executable,
)


class PuckeringDVRLauncher(QObject):
    """
    Asynchronous launcher for the vendored path-DVR workflow.
    """

    finished = Signal(bool, str)

    def __init__(self, workdir: Path, repo_root: Path, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.repo_root = Path(repo_root)
        self.dvr_root = self.repo_root / "puckering_dvr"
        self.script = self.dvr_root / "scripts" / "mw_path_dvr.py"
        self.fortran_bridge = self.dvr_root / "scripts" / "fortran_bridge" / "run_fortran_dvr.py"
        try:
            self.fortran_exe = resolve_dvr_executable(self.repo_root)
        except FileNotFoundError:
            self.fortran_exe = self.repo_root / "bin" / "path_dvr.x"
        self.process = QProcess(self)

    def start_path_analysis(
        self,
        log_path: Path,
        outdir: Path,
        figdir: Path,
        prefix: str,
        boundary: str,
        solver: str,
        compute_rotconst: bool = True,
        label_cremer_pople: bool = True,
    ):
        log_path = Path(log_path)
        if not log_path.exists():
            self.finished.emit(False, f"Gaussian log not found: {log_path}")
            return
        if not self.script.exists():
            self.finished.emit(False, f"Path DVR backend not found: {self.script}")
            return

        request = DVRRequest(
            repo_root=self.repo_root,
            log_path=log_path,
            outdir=Path(outdir),
            figdir=Path(figdir),
            prefix=prefix,
            boundary=boundary,
            solver=solver,
            compute_rotconst=compute_rotconst,
            label_cremer_pople=label_cremer_pople,
        )
        args = build_path_analysis_args(request)
        if is_fortran_solver(solver):
            self._start_fortran_workflow(request, args)
        else:
            self._start_process(sys.executable, args)

    def build_path_analysis_args(
        self,
        log_path: Path,
        outdir: Path,
        figdir: Path,
        prefix: str,
        boundary: str,
        solver: str,
        compute_rotconst: bool = True,
        label_cremer_pople: bool = True,
        check_only: bool = False,
    ) -> list[str]:
        return build_path_analysis_args(
            DVRRequest(
                repo_root=self.repo_root,
                log_path=Path(log_path),
                outdir=Path(outdir),
                figdir=Path(figdir),
                prefix=prefix,
                boundary=boundary,
                solver=solver,
                compute_rotconst=compute_rotconst,
                label_cremer_pople=label_cremer_pople,
                check_only=check_only,
            )
        )

    def _start_fortran_workflow(
        self,
        request: DVRRequest,
        python_args: list[str],
    ) -> None:
        if not self.fortran_bridge.exists():
            self.finished.emit(False, f"Fortran DVR bridge not found: {self.fortran_bridge}")
            return
        if not self.fortran_exe.exists():
            self.finished.emit(False, f"Fortran DVR executable not found: {self.fortran_exe}")
            return

        bridge_args = build_fortran_bridge_args(request, self.fortran_exe)
        shell = build_fortran_shell_command(request, python_args, bridge_args)
        self._start_process("/bin/sh", ["-lc", shell])

    def _start_process(self, executable: str, args: list[str]):
        self.process.setWorkingDirectory(str(self.dvr_root))
        env = QProcess.systemEnvironment()
        py_path = str(self.dvr_root)
        env = [e for e in env if not e.startswith("PYTHONPATH=")]
        env.append("PYTHONPATH=" + py_path)
        self.process.setEnvironment(env)

        self.process.finished.connect(self._on_finished)
        self.process.errorOccurred.connect(self._on_error)
        self.process.start(executable, args)

        if not self.process.waitForStarted(3000):
            self.finished.emit(False, "Failed to start puckering DVR process")

    def _on_finished(self, exit_code, exit_status):
        stdout = bytes(self.process.readAllStandardOutput()).decode(errors="replace").strip()
        stderr = bytes(self.process.readAllStandardError()).decode(errors="replace").strip()
        detail = "\n".join(part for part in (stdout, stderr) if part)
        if exit_code == 0:
            msg = "Path DVR completed"
            if detail:
                msg += "\n\n" + detail
            self.finished.emit(True, msg)
        else:
            msg = f"Path DVR failed (exit_code={exit_code})"
            if detail:
                msg += "\n\n" + detail
            self.finished.emit(False, msg)

    def _on_error(self, error):
        self.finished.emit(False, f"Path DVR process error: {error}")
