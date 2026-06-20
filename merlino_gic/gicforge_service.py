from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
from typing import Iterable

from merlino_core import build_run_manifest, sha256_file, write_manifest
from merlino_fortran import resolve_backend


GICFORGE_OUTPUTS = ("gicforge.out", "provout", "gauin", "msrin", "VPT2in", "bmat.out")


class GICForgeError(RuntimeError):
    """Raised when the GICForge backend cannot be executed successfully."""


@dataclass(frozen=True)
class GICForgeResult:
    workdir: Path
    executable: Path
    logfile: Path
    files: dict[str, Path]
    manifest: Path


def run_gicforge(
    workdir: Path,
    *,
    executable: Path | None = None,
    output_names: Iterable[str] = GICFORGE_OUTPUTS,
) -> GICForgeResult:
    """Run GICForge in `workdir` and write a normalized manifest."""
    run_dir = Path(workdir)
    run_dir.mkdir(parents=True, exist_ok=True)
    exe = Path(executable) if executable is not None else resolve_backend("gicforge")
    logfile = run_dir / "gicforge.log"

    try:
        with logfile.open("w", encoding="utf-8") as log:
            subprocess.run(
                [str(exe)],
                cwd=run_dir,
                stdout=log,
                stderr=log,
                check=True,
            )
    except FileNotFoundError as exc:
        raise GICForgeError(f"Executable not found: {exe}") from exc
    except subprocess.CalledProcessError as exc:
        raise GICForgeError(f"GICForge failed, see {logfile}") from exc

    _copy_legacy_report(run_dir)
    files = _collect_outputs(run_dir, output_names)
    manifest = _write_gicforge_manifest(run_dir, exe, logfile, files)
    return GICForgeResult(
        workdir=run_dir,
        executable=exe,
        logfile=logfile,
        files=files,
        manifest=manifest,
    )


def _copy_legacy_report(run_dir: Path) -> None:
    legacy_output = run_dir / "provout"
    readable_output = run_dir / "gicforge.out"
    if legacy_output.exists():
        readable_output.write_text(
            legacy_output.read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8",
        )


def _collect_outputs(run_dir: Path, output_names: Iterable[str]) -> dict[str, Path]:
    files: dict[str, Path] = {}
    for name in output_names:
        path = run_dir / name
        if path.exists():
            files[name] = path
    return files


def _optional_checksum(path: Path) -> str | None:
    return sha256_file(path) if path.exists() else None


def _write_gicforge_manifest(
    run_dir: Path,
    executable: Path,
    logfile: Path,
    files: dict[str, Path],
) -> Path:
    input_checksums = {
        name: checksum
        for name in ("provin", "xyzin")
        if (checksum := _optional_checksum(run_dir / name)) is not None
    }
    output_checksums = {
        name: sha256_file(path)
        for name, path in sorted(files.items())
        if path.is_file()
    }
    legacy_manifest = {
        "workflow": "gicforge",
        "executable": str(executable),
        "logfile": str(logfile),
        "inputs": input_checksums,
        "outputs": {name: str(path) for name, path in sorted(files.items())},
        "output_sha256": output_checksums,
    }
    manifest = build_run_manifest(
        workflow="gicforge",
        status="completed",
        run_dir=run_dir,
        inputs={name: run_dir / name for name in ("provin", "xyzin") if (run_dir / name).exists()},
        outputs=files,
        backend={"name": "gicforge", "executable": str(executable), "logfile": str(logfile)},
    ).to_dict()
    manifest.update({"legacy": legacy_manifest})
    return write_manifest(run_dir / "gicforge_manifest.json", manifest)
