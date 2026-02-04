"""
Z-matrix reader for Merlino 3.0.

Behavior:
- Reads a Z-matrix file
- Converts it to Cartesian coordinates using zmat_to_xyz
- Replaces ONLY the XYZ block in working/xyzin
- Preserves all existing sections
- Does not return anything
"""

from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Union

from .xyzin_utils import replace_xyz_block
from .zmat_to_xyz import zmat_to_xyz


# ============================================================
# Z-matrix data structures (unchanged)
# ============================================================

@dataclass
class Param:
    name: str
    value: float | None = None


@dataclass
class AngleSpec:
    kind: int
    value: Union[float, Param]


@dataclass
class ZMatrix:
    atoms: List[str]
    refs: List[List[int]]
    bond_length: List[Union[float, Param]]
    bond_angle:  List[Union[float, Param]]
    dihedral:    List[AngleSpec]
    variables: Dict[str, Param] = field(default_factory=dict)

    def resolve(self):
        def _val(x):
            if isinstance(x, Param):
                if x.value is None:
                    raise ValueError(f"Unresolved Z-matrix parameter: {x.name}")
                return x.value
            return x

        return ZMatrix(
            atoms=self.atoms,
            refs=self.refs,
            bond_length=[_val(x) for x in self.bond_length],
            bond_angle=[_val(x) for x in self.bond_angle],
            dihedral=[AngleSpec(d.kind, _val(d.value)) for d in self.dihedral],
            variables=self.variables,
        )


# ============================================================
# Utilities
# ============================================================

def _parse_value(token: str, variables: Dict[str, Param]):
    try:
        return float(token)
    except ValueError:
        if token not in variables:
            variables[token] = Param(token)
        return variables[token]


# ============================================================
# Reader (Merlino 3.0 compliant)
# ============================================================

def read_zmat(path):
    """
    Read a Z-matrix file and replace the XYZ block in xyzin.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    lines = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip()
    ]

    geom_lines = []
    param_lines = []

    # --------------------------------------------------
    # Split geometry and parameters
    # --------------------------------------------------
    for line in lines:
        # Explicit PARAM section
        if line.lower().startswith("param"):
            continue

        # Implicit parameter definition (no PARAM keyword)
        if "=" in line:
            param_lines.append(line)
            continue

        # Geometry line
        geom_lines.append(line)

    # --------------------------------------------------
    # Parse parameters
    # --------------------------------------------------
    variables: Dict[str, Param] = {}
    for line in param_lines:
        if "=" in line:
            k, v = line.split("=", 1)
            if k.strip() not in variables:
                variables[k.strip()] = Param(k.strip())
            variables[k.strip()].value = float(v)

    atoms, refs, bond_length, bond_angle, dihedral = [], [], [], [], []

    # --------------------------------------------------
    # Parse geometry
    # --------------------------------------------------
    for i, line in enumerate(geom_lines):
        f = line.split()
        atoms.append(f[0])

        if i == 0:
            refs.append([0, 0, 0, 0])
            bond_length.append(0.0)
            bond_angle.append(0.0)
            dihedral.append(AngleSpec(0, 0.0))

        elif i == 1:
            refs.append([0, int(f[1]), 0, 0])
            bond_length.append(_parse_value(f[2], variables))
            bond_angle.append(0.0)
            dihedral.append(AngleSpec(0, 0.0))

        elif i == 2:
            refs.append([0, int(f[1]), int(f[3]), 0])
            bond_length.append(_parse_value(f[2], variables))
            bond_angle.append(_parse_value(f[4], variables))
            dihedral.append(AngleSpec(0, 0.0))

        else:
            refs.append([0, int(f[1]), int(f[3]), int(f[5])])
            bond_length.append(_parse_value(f[2], variables))
            bond_angle.append(_parse_value(f[4], variables))
            dihedral.append(AngleSpec(0, _parse_value(f[6], variables)))

    zmat = ZMatrix(
        atoms=atoms,
        refs=refs,
        bond_length=bond_length,
        bond_angle=bond_angle,
        dihedral=dihedral,
        variables=variables,
    ).resolve()

    coords = zmat_to_xyz(zmat)

    symbols = [a for a in zmat.atoms if a not in ("XX", "X", "-1")]

    xyz_block = [
        str(len(symbols)),
        f"Z-matrix input: {path.name}",
    ]

    for s, (x, y, z) in zip(symbols, coords):
        xyz_block.append(f"{s:2s} {x:15.8f} {y:15.8f} {z:15.8f}")

    replace_xyz_block(xyz_block)
