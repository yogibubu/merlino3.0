from .base_launcher import BaseLauncher, LauncherResult
from merlino_gic import GICForgeError, run_gicforge


class GICForgeLauncher(BaseLauncher):
    executable = "gicforge.x"

    def run(self) -> LauncherResult:
        try:
            result = run_gicforge(self.workdir)
        except GICForgeError as exc:
            return LauncherResult(False, {}, str(exc))

        return LauncherResult(
            True,
            result.files | {"gicforge_manifest.json": result.manifest},
            "GICForge completed successfully",
        )
