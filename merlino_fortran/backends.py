from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

from merlino_core.paths import repo_root


@dataclass(frozen=True)
class BackendSpec:
    """Description of a compiled Fortran backend."""

    name: str
    executable: str
    source_dir: str
    build_script: str = "compile_MAC"
    aliases: tuple[str, ...] = ()

    def candidates(self, root: Path) -> list[Path]:
        names = (self.executable, *self.aliases)
        return [root / "bin" / name for name in names]

    def build_command(self, root: Path) -> list[str]:
        return [str(root / self.source_dir / self.build_script)]


BACKENDS: dict[str, BackendSpec] = {
    "gicforge": BackendSpec(
        name="gicforge",
        executable="gicforge.x",
        aliases=("prova.x",),
        source_dir="fortran/gicforge",
    ),
    "dvr": BackendSpec(
        name="dvr",
        executable="path_dvr.x",
        source_dir="fortran/dvr",
    ),
}


def backend_executable(name: str, root: Path | None = None) -> str:
    if name not in BACKENDS:
        known = ", ".join(sorted(BACKENDS))
        raise KeyError(f"Unknown Fortran backend {name!r}; known backends: {known}")
    return BACKENDS[name].executable


def resolve_backend(name: str, root: Path | None = None) -> Path:
    """Resolve a backend executable from `bin/` first, then from PATH."""
    if name not in BACKENDS:
        known = ", ".join(sorted(BACKENDS))
        raise KeyError(f"Unknown Fortran backend {name!r}; known backends: {known}")
    spec = BACKENDS[name]
    root = Path(root) if root is not None else repo_root()

    for candidate in spec.candidates(root):
        if candidate.exists():
            return candidate

    for exe_name in (spec.executable, *spec.aliases):
        found = shutil.which(exe_name)
        if found:
            return Path(found)

    candidates = ", ".join(str(path) for path in spec.candidates(root))
    raise FileNotFoundError(
        f"Fortran backend {name!r} not found. Tried {candidates} and PATH."
    )
