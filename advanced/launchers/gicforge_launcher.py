import subprocess
from .base_launcher import BaseLauncher, LauncherResult
from merlino_fortran import resolve_backend


class GICForgeLauncher(BaseLauncher):
    executable = "gicforge.x"

    def run(self) -> LauncherResult:
        # self.workdir IS the working directory for GICForge.
        run_dir = self.workdir

        logfile = run_dir / "gicforge.log"

        try:
            try:
                exec_path = resolve_backend("gicforge")
            except FileNotFoundError:
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
                False, {}, f"GICForge failed, see {logfile}"
            )

        legacy_output = run_dir / "provout"
        readable_output = run_dir / "gicforge.out"
        if legacy_output.exists():
            readable_output.write_text(
                legacy_output.read_text(encoding="utf-8", errors="replace"),
                encoding="utf-8",
            )

        files = {}
        for name in ["gicforge.out", "provout", "gauin", "msrin", "VPT2in", "bmat.out"]:
            path = run_dir / name
            if path.exists():
                files[name] = path

        return LauncherResult(
            True,
            files,
            "GICForge completed successfully",
        )
