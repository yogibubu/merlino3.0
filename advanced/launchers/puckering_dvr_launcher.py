from __future__ import annotations

from pathlib import Path
import sys

from PySide6.QtCore import QObject, QProcess, Signal


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

        prefix = prefix or "puckering_dvr"
        args = self.build_path_analysis_args(
            log_path,
            outdir,
            figdir,
            prefix,
            boundary,
            solver,
            compute_rotconst=compute_rotconst,
            label_cremer_pople=label_cremer_pople,
        )
        if solver in {"fortran-sinc-dvr", "fortran-gaussian"}:
            self._start_fortran_workflow(args, outdir, prefix, boundary, solver)
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
        effective_solver = "sinc-dvr" if solver in {"fortran-sinc-dvr", "fortran-gaussian"} else solver
        args = [
            str(self.script),
            "--gaussian-log",
            str(log_path),
            "--log-selection",
            "last-per-link",
            "--boundary",
            boundary,
            "--solver",
            effective_solver,
            "--outdir",
            str(outdir),
            "--figdir",
            str(figdir),
            "--prefix",
            prefix or "puckering_dvr",
        ]
        if compute_rotconst:
            args.append("--compute-rotconst")
        if label_cremer_pople:
            args.append("--label-cremer-pople")
        if check_only:
            args.append("--check-only")
        return args

    def _start_fortran_workflow(
        self,
        python_args: list[str],
        outdir: Path,
        prefix: str,
        boundary: str,
        solver: str,
    ) -> None:
        if not self.fortran_bridge.exists():
            self.finished.emit(False, f"Fortran DVR bridge not found: {self.fortran_bridge}")
            return
        if not self.fortran_exe.exists():
            self.finished.emit(False, f"Fortran DVR executable not found: {self.fortran_exe}")
            return

        shell = (
            f"{self._quote(sys.executable)} {' '.join(self._quote(arg) for arg in python_args)}"
            f" && {self._quote(sys.executable)} {self._quote(str(self.fortran_bridge))}"
            f" --grid-csv {self._quote(str(Path(outdir) / f'{prefix}_grid.csv'))}"
            f" --exe {self._quote(str(self.fortran_exe))}"
            f" --outdir {self._quote(str(outdir))}"
            f" --prefix {self._quote(prefix)}"
            f" --boundary {self._quote(boundary)}"
        )
        if solver == "fortran-gaussian":
            shell += " --mode gaussian"
        self._start_process("/bin/sh", ["-lc", shell])

    @staticmethod
    def _quote(value: str) -> str:
        return "'" + value.replace("'", "'\"'\"'") + "'"

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
