"""
Molpro output reader (xyzin-only).

FINAL CONTRACT (Merlino 3.0):
- Reads Molpro output
- Writes XYZ block + #BASIC section to working/xyzin
- Parsing preserved
"""

from pathlib import Path
from typing import List
import re

from .xyzin_utils import replace_xyz_block, append_section


# ==================================================
# Geometry parsing (UNCHANGED)
# ==================================================
def _parse_atomic_coordinates(lines: List[str]):
    indices = [
        i for i, l in enumerate(lines)
        if "ATOMIC COORDINATES" in l.upper()
    ]

    if not indices:
        raise RuntimeError("No 'ATOMIC COORDINATES' section found.")

    i = indices[-1] + 1
    n = len(lines)

    while i < n and not lines[i].strip():
        i += 1

    while i < n:
        parts = lines[i].split()
        try:
            float(parts[-1])
            break
        except Exception:
            i += 1

    atoms, coords = [], []

    while i < n:
        line = lines[i].strip()
        if not line:
            break

        parts = line.split()
        if len(parts) < 4:
            break

        if parts[0].isdigit():
            raw_atom = parts[1]
            xyz = parts[2:5]
        else:
            raw_atom = parts[0]
            xyz = parts[1:4]

        m = re.match(r"[A-Za-z]+", raw_atom)
        if not m:
            raise ValueError(f"Invalid atom label: {raw_atom}")

        atom = m.group()
        x, y, z = map(float, xyz)

        atoms.append(atom)
        coords.append((x, y, z))
        i += 1

    if not atoms:
        raise RuntimeError("Failed to parse Molpro coordinates.")

    return atoms, coords


# ==================================================
# Charge / multiplicity parsing (NEW, MINIMAL)
# ==================================================
def _parse_charge_mult(lines: List[str]):
    charge = 0
    multiplicity = 1

    for line in lines:
        # Explicit charge
        if "CHARGE" in line.upper():
            parts = line.replace("=", " ").split()
            try:
                charge = int(parts[parts.index("CHARGE") + 1])
            except Exception:
                pass

        # Spin quantum number
        if "SPIN QUANTUM NUMBER" in line.upper():
            try:
                s = float(line.split("=")[-1])
                multiplicity = int(2 * s + 1)
            except Exception:
                pass

    return charge, multiplicity


# ==================================================
# Public API
# ==================================================
def read_molpro(path):
    """
    Read Molpro output file and update xyzin (XYZ + #BASIC).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    lines = path.read_text(errors="ignore").splitlines()

    atoms, coords = _parse_atomic_coordinates(lines)
    charge, multiplicity = _parse_charge_mult(lines)

    _write_xyzin(atoms, coords, charge, multiplicity, path.name)


# ==================================================
# WRITE XYzin (MODIFIED – MINIMAL)
# ==================================================
def _write_xyzin(atoms, coords, charge, multiplicity, comment):
    """
    Write XYZ block and #BASIC section to xyzin (Merlino 3.0).
    """
    xyz_lines = [
        str(len(atoms)),
        comment,
    ]

    for a, (x, y, z) in zip(atoms, coords):
        xyz_lines.append(
            f"{a:2s} {x:15.8f} {y:15.8f} {z:15.8f}"
        )

    # --- XYZ
    replace_xyz_block(xyz_lines)

    # --- BASIC (override defaults)
    basic_lines = [
        f"charge {charge}",
        f"multiplicity {multiplicity}",
        "group c1",
        "T_K = 298.15",
        "P_ATM = 1.000000",
    ]
    append_section("BASIC", basic_lines)
