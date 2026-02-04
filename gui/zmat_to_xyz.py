import math
from typing import List, Tuple

# ============================================================
# Vector utilities
# ============================================================

def _norm(v):
    return math.sqrt(sum(x * x for x in v))


def _normalize(v):
    n = _norm(v)
    if n == 0.0:
        raise ValueError("Zero-length vector")
    return [x / n for x in v]


def _cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def _add(a, b):
    return [a[i] + b[i] for i in range(3)]


def _sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def _scale(v, s):
    return [s * x for x in v]


# ============================================================
# Z-matrix → XYZ
# ============================================================

def zmat_to_xyz(zmat) -> List[Tuple[float, float, float]]:
    """
    Convert a RESOLVED ZMatrix (numerical) to Cartesian coordinates.

    Supported kinds:
    - kind = 0   : valence angle + dihedral
    - kind = ±1  : valence angle + second valence angle
    - kind = 2   : out-of-plane angle

    Dummy atoms (XX or Z = -1) are removed from the output.
    """

    # Import locale per evitare import circolare

    from .zmat_reader import ZMatrix

    if not isinstance(zmat, ZMatrix):
        raise TypeError("zmat_to_xyz expects a ZMatrix instance")

    coords: List[List[float]] = []

    # --------------------------------------------------------
    # Build coordinates for ALL atoms (including dummies)
    # --------------------------------------------------------
    for i, (iz, r, a, d) in enumerate(
        zip(zmat.refs, zmat.bond_length, zmat.bond_angle, zmat.dihedral)
    ):

        # Atom 1: origin
        if i == 0:
            coords.append([0.0, 0.0, 0.0])
            continue

        # Atom 2: on Z axis
        if i == 1:
            ref = iz[1] - 1
            coords.append([
                coords[ref][0],
                coords[ref][1],
                coords[ref][2] + r
            ])
            continue

        # Atom 3: define a plane
        if i == 2:
            i1 = iz[1] - 1
            i2 = iz[2] - 1

            theta = math.radians(a)

            p1 = coords[i1]
            p2 = coords[i2]

            ez = _normalize(_sub(p1, p2))

            # deterministic perpendicular
            refv = [0.0, 0.0, 1.0] if abs(ez[2]) < 0.9 else [0.0, 1.0, 0.0]
            ex = _normalize(_cross(refv, ez))
            ey = _cross(ez, ex)

            new = _add(
                p1,
                _add(
                    _scale(ez, -r * math.cos(theta)),
                    _scale(ey, r * math.sin(theta))
                )
            )
            coords.append(new)
            continue

        # ----------------------------------------------------
        # General case (i >= 3)
        # ----------------------------------------------------
        i1 = iz[1] - 1
        i2 = iz[2] - 1
        i3 = iz[3] - 1

        p1 = coords[i1]
        p2 = coords[i2]
        p3 = coords[i3]

        # Primary axis (bond direction)
        ez = _normalize(_sub(p1, p2))

        # Secondary reference (physical if possible)
        vref = _sub(p3, p2)
        if _norm(_cross(ez, vref)) > 1e-8:
            ex = _normalize(_cross(vref, ez))
        else:
            # fallback deterministic axis
            refv = [0.0, 0.0, 1.0] if abs(ez[2]) < 0.9 else [0.0, 1.0, 0.0]
            ex = _normalize(_cross(refv, ez))

        ey = _cross(ez, ex)

        theta = math.radians(a)

        # ----------------------------
        # kind = 0 : angle + dihedral
        # ----------------------------
        if d.kind == 0:
            phi = math.radians(d.value)

            new = _add(
                p1,
                _add(
                    _scale(ez, -r * math.cos(theta)),
                    _add(
                        _scale(ex, r * math.sin(theta) * math.cos(phi)),
                        _scale(ey, r * math.sin(theta) * math.sin(phi))
                    )
                )
            )

        # -----------------------------------------
        # kind = ±1 : angle + second valence angle
        # -----------------------------------------
        elif d.kind in (+1, -1):
            theta2 = math.radians(d.value)

            if d.kind == +1:
                ep, eq = ex, ey
            else:
                ep, eq = ey, ex

            new = _add(
                p1,
                _add(
                    _scale(ez, -r * math.cos(theta)),
                    _add(
                        _scale(ep, r * math.sin(theta) * math.cos(theta2)),
                        _scale(eq, r * math.sin(theta) * math.sin(theta2))
                    )
                )
            )

        # -----------------------------------------
        # kind = 2 : out-of-plane angle
        # -----------------------------------------
        elif d.kind == 2:
            oop = math.radians(d.value)

            rin = r * math.cos(oop)
            rout = r * math.sin(oop)

            new = _add(
                p1,
                _add(
                    _scale(ez, -rin * math.cos(theta)),
                    _add(
                        _scale(ex, rin * math.sin(theta)),
                        _scale(ey, rout)
                    )
                )
            )

        else:
            raise ValueError(f"Unsupported dihedral kind: {d.kind}")

        coords.append(new)

    # --------------------------------------------------------
    # Remove dummy atoms
    # --------------------------------------------------------
    clean_coords: List[Tuple[float, float, float]] = []

    for atom, xyz in zip(zmat.atoms, coords):
        if atom == "XX" or atom == -1:
            continue
        clean_coords.append(tuple(xyz))

    return clean_coords


# ============================================================
# Write XYZ
# ============================================================

def write_xyz(zmat, filename: str):
    """
    Write an XYZ file from a RESOLVED ZMatrix.
    """

    from .zmatrix import ZMatrix
    if not isinstance(zmat, ZMatrix):
        raise TypeError("write_xyz expects a ZMatrix instance")

    xyz = zmat_to_xyz(zmat)

    atoms = []
    for atom in zmat.atoms:
        if atom == "XX" or atom == -1:
            continue
        atoms.append(atom)

    if len(atoms) != len(xyz):
        raise RuntimeError("Atom/coordinate mismatch in write_xyz")

    with open(filename, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write("\n")
        for atom, (x, y, z) in zip(atoms, xyz):
            f.write(f"{atom:>2s}  {x:15.8f}  {y:15.8f}  {z:15.8f}\n")

