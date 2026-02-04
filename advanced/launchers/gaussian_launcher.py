from pathlib import Path
from PySide6.QtCore import QObject, QProcess, Signal


class GaussianLauncher(QObject):
    """
    Asynchronous Gaussian launcher using QProcess.
    """

    finished = Signal(bool, str)  # success, message

    def __init__(self, workdir: Path, parent=None):
        super().__init__(parent)
        self.workdir = Path(workdir)
        self.process = QProcess(self)

    def _select_log_path(self):
        """
        Select the most likely active Gaussian log.
        Prefer the most recently modified existing file between gauin.log/gauout.log.
        """
        candidates = []
        for name in ("gauin.log", "gauout.log"):
            p = self.workdir / name
            if p.exists():
                try:
                    st = p.stat()
                except Exception:
                    continue
                candidates.append((st.st_mtime, st.st_size, p))
        if not candidates:
            return self.workdir / "gauin.log"
        candidates.sort(key=lambda t: (t[0], t[1]), reverse=True)
        return candidates[0][2]

    def start(self):
        gauin_gjf = self.workdir / "gauin.gjf"
        gauin_raw = self.workdir / "gauin"
        if not gauin_gjf.exists() and gauin_raw.exists():
            try:
                gauin_gjf.write_text(gauin_raw.read_text())
            except Exception:
                pass

        gauin = gauin_gjf if gauin_gjf.exists() else gauin_raw
        if not gauin.exists():
            self.finished.emit(False, "gauin.gjf or gauin not found")
            return

        # Gaussian executable (development version)
        executable = "gdv"

        self.process.setWorkingDirectory(str(self.workdir))

        self.process.finished.connect(self._on_finished)
        self.process.errorOccurred.connect(self._on_error)

        # Start Gaussian in background
        self.process.start(executable, [str(gauin)])

        if not self.process.waitForStarted(3000):
            self.finished.emit(False, "Failed to start Gaussian process")

    def _on_finished(self, exit_code, exit_status):
        log_path = self._select_log_path()
        if log_path.exists():
            try:
                tail = log_path.read_text(encoding="utf-8", errors="replace")
                ok = "Normal termination of Gaussian" in tail
            except Exception:
                ok = False
            if ok:
                self.finished.emit(True, "Gaussian completed successfully")
                return

        if exit_code == 0 and log_path.exists():
            self.finished.emit(True, f"Gaussian completed (check {log_path.name})")
        else:
            self.finished.emit(
                False,
                f"Gaussian finished with errors (exit_code={exit_code}; see {log_path.name})"
            )

    def _on_error(self, error):
        self.finished.emit(False, f"Gaussian process error: {error}")
