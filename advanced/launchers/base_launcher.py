from pathlib import Path
from typing import Dict, Optional
import shutil


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
    1. merlino3.0/bin/<executable>
    2. system PATH
    """

    executable: str = ""

    def __init__(self, workdir: Path):
        self.workdir = Path(workdir)

    # ------------------------------------------------------------------

    def resolve_executable(self) -> Path:
        """
        Resolve the executable path.

        First look for Merlino internal executables in merlino3.0/bin,
        then fallback to PATH.
        """
        # File position:
        # merlino3.0/advanced/launchers/base_launcher.py
        # parents[0] -> launchers
        # parents[1] -> advanced
        # parents[2] -> merlino3.0
        repo_root = Path(__file__).resolve().parents[2]

        # 1. internal Merlino executable
        local_exec = repo_root / "bin" / self.executable
        if local_exec.exists():
            return local_exec

        # 2. system PATH
        path_exec = shutil.which(self.executable)
        if path_exec:
            return Path(path_exec)

        raise FileNotFoundError(
            f"Executable '{self.executable}' not found in merlino3.0/bin or PATH"
        )

    # ------------------------------------------------------------------

    def run(self) -> LauncherResult:
        raise NotImplementedError

