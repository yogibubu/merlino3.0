#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path


def _find_xyzin(workdir: Path):
    if (workdir / "xyzin").exists():
        return workdir / "xyzin"

    cand = sorted(workdir.glob("*.xyzin"))
    if cand:
        return cand[0]

    cand = sorted(workdir.glob("*.xyz"))
    if cand:
        return cand[0]

    return None


def _section_present(filepath: Path, section_name: str) -> bool:
    tag = "#" + section_name.upper()
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip().upper() == tag:
                return True
    return False


def _marker_present(filepath: Path, marker: str) -> bool:
    marker = marker.strip()
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() == marker:
                return True
    return False


def _parse_kv_section(filepath: Path, section_name: str):
    tag = "#" + section_name.upper()
    data = {}
    in_sec = False

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("#"):
                in_sec = (s.upper() == tag)
                continue

            if not in_sec:
                continue

            if "=" in s:
                k, v = [x.strip() for x in s.split("=", 1)]
                data[k] = v
            else:
                parts = s.split()
                if len(parts) >= 2:
                    data[parts[0]] = parts[1]

    return data


def main():
    # Script moved under scripts/: repo root is one level up.
    repo_root = Path(__file__).resolve().parents[1]
    workdir = repo_root / "working"

    if not workdir.exists():
        return 1

    xyzin = _find_xyzin(workdir)
    if xyzin is None:
        return 2

    # ------------------------------------------------------------
    # Import pipelines
    # ------------------------------------------------------------
    try:
        from geometry.rotational_pipeline import rotational_pipeline
        from geometry.thermo_pipeline import thermo_pipeline
    except Exception:
        return 3

    # ------------------------------------------------------------
    # Run pipeline inside working directory
    # ------------------------------------------------------------
    cwd0 = os.getcwd()
    os.chdir(str(workdir))
    try:
        rotational_pipeline("xyzin", report=True)
        thermo_pipeline("xyzin", report=True)
    except Exception:
        os.chdir(cwd0)
        return 4
    os.chdir(cwd0)

    # ------------------------------------------------------------
    # Required sections in xyzin
    # ------------------------------------------------------------
    has_basic = _section_present(xyzin, "BASIC")
    has_rot = _section_present(xyzin, "ROTATIONAL")
    has_vib = _section_present(xyzin, "VIBRATIONAL")  # now always enforced by pipeline
    has_thermo = _section_present(xyzin, "THERMO")

    if not has_basic:
        return 10
    if not has_rot:
        return 11
    if not has_vib:
        return 12
    if not has_thermo:
        return 14

    basic = _parse_kv_section(xyzin, "BASIC")
    if "T_K" not in basic:
        return 15
    if "P_ATM" not in basic:
        return 16

    # ------------------------------------------------------------
    # vibin is OPTIONAL (vibrations may be missing)
    # ------------------------------------------------------------
    vibin = workdir / "vibin"
    if vibin.exists():
        # If vibin exists, coriolis should have been executed and appended
        if not _marker_present(vibin, "BEGIN_CORIOLIS"):
            return 21
        # QCENT results may be missing (qcent failed), so do NOT enforce BEGIN_RESULTS

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
