import subprocess
from pathlib import Path
from .base_launcher import BaseLauncher, LauncherResult


class ProvaLauncher(BaseLauncher):
    executable = "prova.x"

    def run(self) -> LauncherResult:
        # self.workdir IS the working directory for prova.x
        run_dir = self.workdir

        logfile = run_dir / "prova.log"

        try:
            exec_path = self.resolve_executable()

            with logfile.open("w") as log:
                subprocess.run(
                    [str(exec_path)],
                    cwd=run_dir,
                    stdout=log,
                    stderr=log,
                    check=True,
                )

        except FileNotFoundError as e:
            return LauncherResult(
                False, {}, f"Executable not found: {e}"
            )

        except subprocess.CalledProcessError:
            return LauncherResult(
                False, {}, f"prova failed, see {logfile}"
            )

        files = {}
        for name in ["provout", "gauin", "msrin", "VPT2in"]:
            path = run_dir / name
            if path.exists():
                files[name] = path

        return LauncherResult(
            True,
            files,
            "prova completed successfully",
        )

