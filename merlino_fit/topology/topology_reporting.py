import numpy as np
from pathlib import Path

try:
    from survibfit.primitives import build_primitives
    from survibfit.symmetry_classifier import group_label as _group_label
    from survibfit.symmetry_detector import (
        is_linear,
        orient_coords,
        symmetry_elements_from_geometry,
    )
    from survibfit.symmetry_global import primitive_permutation
except Exception:  # pragma: no cover - fallback when `survibfit` is not top-level
    from merlino_fit.survibfit.primitives import build_primitives
    from merlino_fit.survibfit.symmetry_classifier import group_label as _group_label
    from merlino_fit.survibfit.symmetry_detector import (
        is_linear,
        orient_coords,
        symmetry_elements_from_geometry,
    )
    from merlino_fit.survibfit.symmetry_global import primitive_permutation


# ============================================================
# Complete periodic table
# ============================================================

_PERIODIC_TABLE = {
    1: "H",   2: "He",
    3: "Li",  4: "Be",  5: "B",   6: "C",   7: "N",   8: "O",
    9: "F",  10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",  16: "S",
    17: "Cl", 18: "Ar",
    19: "K",  20: "Ca", 21: "Sc", 22: "Ti", 23: "V",  24: "Cr",
    25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y",  40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd",
    49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I",  54: "Xe",
    55: "Cs", 56: "Ba",
    57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm",
    63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er",
    69: "Tm", 70: "Yb", 71: "Lu",
    72: "Hf", 73: "Ta", 74: "W",  75: "Re", 76: "Os", 77: "Ir",
    78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn",
    87: "Fr", 88: "Ra",
    89: "Ac", 90: "Th", 91: "Pa", 92: "U",  93: "Np", 94: "Pu",
    95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm",
    101: "Md", 102: "No", 103: "Lr",
    104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs",
    109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn",
    113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts",
    118: "Og",
}


def _element_symbol(Z):
    """Return element symbol from atomic number."""
    return _PERIODIC_TABLE.get(int(Z), f"Z{int(Z)}")


def _uf_find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


def _uf_union(parent, a, b):
    ra = _uf_find(parent, a)
    rb = _uf_find(parent, b)
    if ra != rb:
        parent[rb] = ra


def _primitive_label(p):
    if p.kind == "bond":
        i, j = p.atoms
        return f"({i+1}-{j+1})"
    if p.kind == "angle":
        i, j, k = p.atoms
        return f"({i+1}-{j+1}-{k+1})"
    if p.kind == "dihedral":
        i, j, k, l = p.atoms
        return f"({i+1}-{j+1}-{k+1}-{l+1})"
    return "(" + ",".join(str(a + 1) for a in p.atoms) + ")"


def _symmetry_summary(cg, dg):
    """
    Compute point group and symmetry-equivalent primitive classes.
    Uses the same geometry already available in topology.
    """
    Z = np.asarray(cg.Z, dtype=int)
    coords = np.asarray(cg.coords, dtype=float)
    symbols = [_element_symbol(z) for z in Z]

    coords_oriented = orient_coords(coords, weights=Z.astype(float))
    elements, atom_classes, permutations = symmetry_elements_from_geometry(
        symbols,
        coords_oriented,
        tol=1.0e-3,
        max_n=8,
        auto_max_n=True,
    )
    point_group = _group_label(elements, linear=is_linear(coords_oriented, tol=1.0e-3))

    prims = build_primitives(dg, coords)
    nprim = len(prims)
    parent = list(range(nprim))
    for mapping in permutations:
        perm_idx, _ = primitive_permutation(prims, mapping)
        for i, j in enumerate(perm_idx):
            _uf_union(parent, i, j)

    classes = {}
    for i in range(nprim):
        root = _uf_find(parent, i)
        classes.setdefault(root, []).append(i)

    by_kind = {}
    for idxs in classes.values():
        idxs_sorted = sorted(idxs)
        kind = prims[idxs_sorted[0]].kind
        by_kind.setdefault(kind, []).append(idxs_sorted)

    return {
        "point_group": point_group,
        "atom_classes": [sorted(cls) for cls in atom_classes],
        "primitive_classes_by_kind": by_kind,
        "primitives": prims,
    }


# ============================================================
# Topology reporting
# ============================================================

def print_topology_report(cg, dg, synthons, arom=None, filename="topology.report", ringset=None):
    """
    Write a diagnostic report of the molecular topology to a file.

    Atom indices in the report are 1-based (Merlino convention).
    """

    Z = cg.Z
    coords = cg.coords

    # --- write output into cwd/working ---
    outdir = Path.cwd() / "working"
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / filename

    with open(outfile, "w") as fh:
        charge_source = getattr(synthons, "_charge_source", "Synthons model")
        bo_source = getattr(synthons, "_bond_order_source", "Topology continuous model")

        fh.write("\nTOPOLOGY INPUT SOURCES\n")
        fh.write("----------------------\n")
        fh.write(f"Atomic charges source: {charge_source}\n")
        fh.write(f"Bond orders source: {bo_source}\n")

        # ========================================================
        # Bonded pairs and bond orders
        # ========================================================

        fh.write("\nBOND ORDERS (bonded pairs, heavy atoms only)\n")
        fh.write("------------------------------------------\n")

        for (i, j) in dg.bonds:
            if Z[i] <= 1 or Z[j] <= 1:
                continue

            rij = np.linalg.norm(coords[i] - coords[j])
            zi = _element_symbol(Z[i])
            zj = _element_symbol(Z[j])

            # >>> CORRECT bond order <<<
            bo = synthons.bond_order(i, j)

            fh.write(
                f"({i+1:2d},{j+1:2d})  "
                f"{zi:>2s}–{zj:<2s}   "
                f"r = {rij:6.3f} Å   "
                f"BO = {bo:6.3f}\n"
            )

        # ========================================================
        # Atomic synthons
        # ========================================================

        fh.write("\nATOMIC SYNTHONS (heavy atoms)\n")
        fh.write("----------------------------\n")

        for i in range(len(Z)):
            if Z[i] <= 1:
                continue

            zi = _element_symbol(Z[i])

            cna = synthons.cna(i)
            edom = synthons._electron_domains(i)
            cn_strain = int(round(edom))

            fh.write(f"\nAtom {i+1:2d}  {zi}  (Z={Z[i]})\n")
            fh.write(f"  Zeff              = {synthons.Zeff(i):8.3f}\n")
            fh.write(f"  CNA               = {cna:8.3f}\n")
            fh.write(f"  Electron domains  = {edom:8.3f}\n")
            fh.write(f"  CN (strain)       = {cn_strain:8d}\n")
            fh.write(f"  Charge            = {synthons.charge(i):8.3f}\n")
            fh.write(f"  Covalency         = {synthons.covalency(i):8.3f}\n")
            fh.write(f"  Delocalization    = {synthons.delocalization(i):8.3f}\n")
            fh.write(f"  Strain            = {synthons.strain(i):8.3f}\n")

        # ========================================================
        # Rings (optional)
        # ========================================================
        if ringset is not None:
            fh.write("\nRINGS\n")
            fh.write("-----\n")
            rings = list(getattr(ringset, "rings", []))
            if rings:
                for ring in rings:
                    atoms = [a + 1 for a in getattr(ring, "atoms", [])]
                    fh.write(
                        f"Ring {int(getattr(ring, 'index', 0)) + 1}: "
                        f"size={len(atoms)} atoms={atoms}\n"
                    )
            else:
                fh.write("No rings detected\n")

        # ========================================================
        # Aromaticity (optional)
        # ========================================================

        if arom is not None:
            fh.write("\nAROMATICITY\n")
            fh.write("-----------\n")

            if arom.aromatic_atoms:
                fh.write(
                    "Aromatic atoms: "
                    f"{sorted([i+1 for i in arom.aromatic_atoms])}\n"
                )
            else:
                fh.write("Aromatic atoms: none\n")

            if arom.aromatic_bonds:
                fh.write(
                    "Aromatic bonds: "
                    f"{sorted([(i+1, j+1) for (i, j) in arom.aromatic_bonds])}\n"
                )
            else:
                fh.write("Aromatic bonds: none\n")

        # ========================================================
        # Global symmetry and equivalent parameters
        # ========================================================
        try:
            sym = _symmetry_summary(cg, dg)
            fh.write("\nGLOBAL SYMMETRY\n")
            fh.write("---------------\n")
            fh.write(f"Point group: {sym['point_group']}\n")

            atom_classes = [
                [i + 1 for i in cls] for cls in sym["atom_classes"] if len(cls) > 1
            ]
            if atom_classes:
                fh.write("Equivalent atom classes (1-based):\n")
                for icls, cls in enumerate(atom_classes, start=1):
                    fh.write(f"  class {icls:2d}: {cls}\n")
            else:
                fh.write("Equivalent atom classes: none (C1-like partition)\n")

            fh.write("\nEQUIVALENT INTERNAL-PARAMETER CLASSES\n")
            fh.write("------------------------------------\n")
            prims = sym["primitives"]
            by_kind = sym["primitive_classes_by_kind"]
            for kind in ("bond", "angle", "dihedral", "out_of_plane", "linear_bend"):
                groups = [g for g in by_kind.get(kind, []) if len(g) > 1]
                fh.write(f"{kind}: {len(groups)} equivalent classes\n")
                for ig, g in enumerate(groups, start=1):
                    labels = [_primitive_label(prims[idx]) for idx in g]
                    fh.write(f"  class {ig:2d}: {labels}\n")
        except Exception as exc:
            fh.write("\nGLOBAL SYMMETRY\n")
            fh.write("---------------\n")
            fh.write(f"Symmetry analysis unavailable: {exc}\n")
