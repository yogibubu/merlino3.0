"""
MRCC output reader (xyzin-only).

FINAL CONTRACT (Merlino 3.0):
- Reads MRCC output
- Writes XYZ block + #BASIC section to working/xyzin
- Parsing preserved
"""

from pathlib import Path
from typing import List

from .xyzin_utils import replace_xyz_block, append_section


# ==================================================
# Geometry parsing (UNCHANGED)
# ==================================================
def _parse_cartesian_coordinates(lines: List[str]):
    indices = [
        i for i, l in enumerate(lines)
        if "CARTESIAN COORDINATES" in l.upper()
    ]

    if not indices:
        raise RuntimeError("No Cartesian coordinates section found.")

    i = indices[-1] + 2
    n = len(lines)

    atoms, coords = [], []

    while i < n:
        line = lines[i].strip()
        if not line:
            break

        parts = line.split()
        if len(parts) < 4:
            break

        atom = parts[0]
        x, y, z = map(float, parts[1:4])

        atoms.append(atom)
        coords.append((x, y, z))
        i += 1

    if not atoms:
        raise RuntimeError("Failed to parse MRCC coordinates.")

    return atoms, coords


# ==================================================
# Charge / multiplicity parsing (NEW, MINIMAL)
# ==================================================
def _parse_charge_mult(lines: List[str]):
    charge = 0
    multiplicity = 1

    for line in lines:
        u = line.upper()

        if "CHARGE OF THE SYSTEM" in u:
            try:
                charge = int(line.split(":")[-1].strip())
            except Exception:
                pass

        if "SPIN MULTIPLICITY" in u:
            try:
                multiplicity = int(line.split(":")[-1].strip())
            except Exception:
                pass

    return charge, multiplicity


# ==================================================
# Public API
# ==================================================
def read_mrcc(path):
    """
    Read MRCC output file and update xyzin (XYZ + #BASIC).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    lines = path.read_text(errors="ignore").splitlines()

    atoms, coords = _parse_cartesian_coordinates(lines)
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
