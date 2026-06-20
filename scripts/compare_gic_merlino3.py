from __future__ import annotations

import argparse
import filecmp
import json
import shutil
import subprocess
import sys
from pathlib import Path


OUTPUTS = ("gauin", "provout", "gicforge.out")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare Merlino3 and Merlino4 GICForge outputs on one fixture.")
    parser.add_argument("--merlino3", type=Path, default=Path("/Users/vincenzobarone/merlino3.0"))
    parser.add_argument("--merlino4", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--fixture", type=Path, default=None, help="Directory containing at least provin")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON report path")
    parser.add_argument("--strict", action="store_true", help="Return failure if baseline/executables/fixture are missing")
    args = parser.parse_args(argv)

    report = compare_gic_outputs(args.merlino3, args.merlino4, args.fixture)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] == "passed":
        return 0
    if args.strict:
        return 1
    return 0


def compare_gic_outputs(merlino3: Path, merlino4: Path, fixture: Path | None) -> dict:
    m3 = Path(merlino3)
    m4 = Path(merlino4)
    exe3 = _gic_executable(m3)
    exe4 = _gic_executable(m4)
    if not m3.exists():
        return _skipped(f"Merlino3 baseline not found: {m3}")
    if exe3 is None:
        return _skipped(f"Merlino3 GIC executable not found below {m3}")
    if exe4 is None:
        return _skipped(f"Merlino4 GIC executable not found below {m4}")
    if fixture is None or not (Path(fixture) / "provin").exists():
        return _skipped("fixture with provin was not provided")

    work_root = m4 / ".tmp_gic_regression"
    if work_root.exists():
        shutil.rmtree(work_root)
    run3 = work_root / "merlino3"
    run4 = work_root / "merlino4"
    shutil.copytree(fixture, run3)
    shutil.copytree(fixture, run4)
    try:
        res3 = _run_gic(exe3, run3)
        res4 = _run_gic(exe4, run4)
        compared = {}
        mismatches = []
        for name in OUTPUTS:
            p3 = run3 / name
            p4 = run4 / name
            if p3.exists() and p4.exists():
                same = _normalized_text(p3) == _normalized_text(p4)
                compared[name] = same
                if not same:
                    mismatches.append(name)
        status = "passed" if res3 == 0 and res4 == 0 and compared and not mismatches else "failed"
        return {
            "status": status,
            "merlino3": str(m3),
            "merlino4": str(m4),
            "fixture": str(fixture),
            "return_codes": {"merlino3": res3, "merlino4": res4},
            "compared": compared,
            "mismatches": mismatches,
        }
    finally:
        shutil.rmtree(work_root, ignore_errors=True)


def _gic_executable(root: Path) -> Path | None:
    for rel in ("bin/gicforge.x", "bin/prova.x", "fortran/gicforge/gicforge.x", "fortran/gicforge/prova.x"):
        path = root / rel
        if path.exists():
            return path
    return None


def _run_gic(exe: Path, cwd: Path) -> int:
    try:
        return subprocess.run([str(exe)], cwd=cwd, timeout=60, stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode
    except Exception:
        return 999


def _normalized_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    return "\n".join(line.rstrip() for line in text.splitlines() if line.strip()) + "\n"


def _skipped(reason: str) -> dict:
    return {"status": "skipped", "reason": reason}


if __name__ == "__main__":
    raise SystemExit(main())
