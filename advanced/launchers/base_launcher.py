from pathlib import Path
from typing import Dict, Optional
import shutil

from merlino_core.paths import repo_root


class LauncherResult:
    def __init__(
        self,
        success: bool,
        files: Dict[str, Path],
        message: str = "",
        process: Optional[object] = None,
    ):
        self.success = success
        self.files = files
        self.message = message
        self.process = process


class BaseLauncher:
    """
    Base class for all Merlino launchers.

    Executable resolution order:
    1. Merlino repository bin/<executable>
    2. system PATH
    """

    executable: str = ""

    def __init__(self, workdir: Path):
        self.workdir = Path(workdir)

    # ------------------------------------------------------------------

    def resolve_executable(self) -> Path:
        """
        Resolve the executable path.

        First look for Merlino internal executables in bin/, then fallback to
        PATH. New Fortran-specific wrappers live in `merlino_fortran`; this
        method remains for generic launcher compatibility.
        """
        root = repo_root(__file__)

        # 1. internal Merlino executable
        local_exec = root / "bin" / self.executable
        if local_exec.exists():
            return local_exec

        # 2. system PATH
        path_exec = shutil.which(self.executable)
        if path_exec:
            return Path(path_exec)

        raise FileNotFoundError(
            f"Executable '{self.executable}' not found in {root / 'bin'} or PATH"
        )

    # ------------------------------------------------------------------

    def run(self) -> LauncherResult:
        raise NotImplementedError
