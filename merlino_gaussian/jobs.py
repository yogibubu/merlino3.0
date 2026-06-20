from __future__ import annotations

from pathlib import Path


GAUSSIAN_EXECUTABLE = "gdv"
NORMAL_TERMINATION_MARKER = "Normal termination of Gaussian"
LOG_CANDIDATES = ("gauin.log", "gauout.log")


class GaussianInputError(RuntimeError):
    """Raised when no runnable Gaussian input can be found."""


def ensure_gjf_input(workdir: Path) -> Path:
    """Return a Gaussian input path, copying `gauin` to `gauin.gjf` if needed."""
    workdir = Path(workdir)
    gauin_gjf = workdir / "gauin.gjf"
    gauin_raw = workdir / "gauin"
    if not gauin_gjf.exists() and gauin_raw.exists():
        gauin_gjf.write_text(
            gauin_raw.read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8",
        )
    if gauin_gjf.exists():
        return gauin_gjf
    if gauin_raw.exists():
        return gauin_raw
    raise GaussianInputError("gauin.gjf or gauin not found")


def select_latest_log(workdir: Path) -> Path:
    """Select the most likely active Gaussian log in a work directory."""
    workdir = Path(workdir)
    candidates: list[tuple[float, int, Path]] = []
    for name in LOG_CANDIDATES:
        path = workdir / name
        if path.exists():
            try:
                stat = path.stat()
            except OSError:
                continue
            candidates.append((stat.st_mtime, stat.st_size, path))
    if not candidates:
        return workdir / "gauin.log"
    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return candidates[0][2]


def gaussian_completed_normally(log_path: Path) -> bool:
    if not Path(log_path).exists():
        return False
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    return NORMAL_TERMINATION_MARKER in text


def gaussian_completion_message(workdir: Path, exit_code: int) -> tuple[bool, str]:
    """Return `(success, message)` for a finished Gaussian process."""
    log_path = select_latest_log(workdir)
    if gaussian_completed_normally(log_path):
        return True, "Gaussian completed successfully"
    if exit_code == 0 and log_path.exists():
        return True, f"Gaussian completed (check {log_path.name})"
    return False, f"Gaussian finished with errors (exit_code={exit_code}; see {log_path.name})"
