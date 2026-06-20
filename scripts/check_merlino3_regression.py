from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Check Merlino4 against local Merlino3 frozen baseline.")
    parser.add_argument("--merlino3", type=Path, default=Path("/Users/vincenzobarone/merlino3.0"))
    parser.add_argument("--merlino4", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--skip-tests", action="store_true")
    args = parser.parse_args(argv)

    failures: list[str] = []
    if not args.merlino3.exists():
        failures.append(f"Merlino3 baseline not found: {args.merlino3}")
    for rel in ("app.py", "advanced/advanced_window.py", "advanced/dvr_window.py", "bin/gicforge.x"):
        if not (args.merlino4 / rel).exists():
            failures.append(f"Merlino4 missing expected entry point: {rel}")
    for rel in ("app.py", "advanced/advanced_window.py"):
        if args.merlino3.exists() and not (args.merlino3 / rel).exists():
            failures.append(f"Merlino3 baseline missing expected file: {rel}")

    if not args.skip_tests:
        result = subprocess.run(["./freeze_check.sh"], cwd=args.merlino4)
        if result.returncode != 0:
            failures.append("Merlino4 freeze_check.sh failed")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    print(f"Merlino3 baseline: {args.merlino3}")
    print(f"Merlino4 checked: {args.merlino4}")
    print("No baseline-entry regression detected.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
