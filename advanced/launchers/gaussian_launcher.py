from pathlib import Path
from PySide6.QtCore import QObject, QProcess, Signal
from merlino_gaussian import (
    GAUSSIAN_EXECUTABLE,
    GaussianInputError,
    ensure_gjf_input,
    gaussian_completion_message,
    select_latest_log,
)


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
        return select_latest_log(self.workdir)

    def start(self):
        try:
            gauin = ensure_gjf_input(self.workdir)
        except GaussianInputError as exc:
            self.finished.emit(False, str(exc))
            return

        self.process.setWorkingDirectory(str(self.workdir))

        self.process.finished.connect(self._on_finished)
        self.process.errorOccurred.connect(self._on_error)

        # Start Gaussian in background
        self.process.start(GAUSSIAN_EXECUTABLE, [str(gauin)])

        if not self.process.waitForStarted(3000):
            self.finished.emit(False, "Failed to start Gaussian process")

    def _on_finished(self, exit_code, exit_status):
        success, message = gaussian_completion_message(self.workdir, int(exit_code))
        self.finished.emit(success, message)

    def _on_error(self, error):
        self.finished.emit(False, f"Gaussian process error: {error}")
