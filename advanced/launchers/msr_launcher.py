# advanced/launchers/msr_launcher.py

import subprocess
from pathlib import Path
from .base_launcher import BaseLauncher, LauncherResult


class MSRLauncher(BaseLauncher):
    executable = "msr.x"

    def run(self) -> LauncherResult:
        logfile = self.workdir / "msr.log"

        try:
            with logfile.open("w") as log:
                subprocess.run(
                    [self.executable],
                    cwd=self.workdir,
                    stdout=log,
                    stderr=log,
                    check=True,
                )
        except subprocess.CalledProcessError:
            return LauncherResult(
                False, {}, f"MSR failed, see {logfile.name}"
            )

        outputs = {}
        out = self.workdir / "msrout"
        if out.exists():
            outputs["msrout"] = out

        return LauncherResult(
            True,
            outputs,
            "MSR completed successfully",
        )

