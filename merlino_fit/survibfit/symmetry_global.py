from __future__ import annotations

import numpy as np
import re

from topology.elements import atomic_symbol

from .symmetry_classifier import (
    group_label as _group_label,
    highest_cn_axis as _highest_cn_axis,
    pick_generators as _pick_generators,
    spherical_top_guess as _spherical_top_guess,
)
from .symmetry_detector import (
    is_linear,
    orient_coords,
    symmetry_elements_from_geometry,
)


def irrep_characters_for_operations(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    """Return real irrep projectors for the detected operation sequence.

    The returned character vectors are scaled by the irrep dimension.  This
    keeps the existing projection formula `(1/h) sum chi(R) R` valid for both
    one-dimensional and degenerate irreps, and makes the vibrational target
    count equal to the dimension of the projected subspace.
    """
    canonical = [_canonical_operation_label(label) for label in labels]
    if canonical == ["E"]:
        return [("A", np.ones(1))]
    group = _group_label([(label, None, 0.0) for label in labels])
    if group == "Cs":
        return _cs_irreps(canonical)
    if group == "Ci":
        return _ci_irreps(canonical)
    if group == "C2":
        return _c2_irreps(canonical)
    if group == "C2v":
        return _c2v_irreps(canonical)
    if group == "D2":
        return _d2_irreps(canonical)
    if group == "C2h":
        return _c2h_irreps(canonical)
    if group == "D2h":
        return _d2h_irreps(canonical)
    general = _cyclic_family_irreps(group, canonical)
    if general:
        return general
    return []


def _canonical_operation_label(label: str) -> str:
    text = str(label)
    if text == "E":
        return "E"
    if text == "i":
        return "i"
    if text.startswith("sigma"):
        return text
    match = re.match(r"C(\d+)([xyz])\^(\d+)", text)
    if match:
        n, axis, power = match.groups()
        if int(n) == 2:
            return f"C2{axis}"
        return f"C{n}{axis}^{power}"
    if text.startswith("C2_xy"):
        return "C2_perp"
    if text.startswith("S"):
        return text
    return text


def _chars(labels: list[str], values: dict[str, float]) -> np.ndarray:
    return np.array([values.get(label, 1.0) for label in labels], dtype=float)


def _cs_irreps(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    return [("A'", _chars(labels, {"sigma_xy": 1.0, "sigma_xz": 1.0, "sigma_yz": 1.0})),
            ("A''", _chars(labels, {"sigma_xy": -1.0, "sigma_xz": -1.0, "sigma_yz": -1.0}))]


def _ci_irreps(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    return [("Ag", _chars(labels, {"i": 1.0})), ("Au", _chars(labels, {"i": -1.0}))]


def _c2_irreps(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    return [("A", _chars(labels, {"C2x": 1.0, "C2y": 1.0, "C2z": 1.0, "C2_perp": 1.0})),
            ("B", _chars(labels, {"C2x": -1.0, "C2y": -1.0, "C2z": -1.0, "C2_perp": -1.0}))]


def _c2v_irreps(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    sigma_labels = [label for label in labels if label.startswith("sigma")]
    preferred = ("sigma_xz", "sigma_yz", "sigma_xy")
    first_sigma = next((label for label in preferred if label in sigma_labels), sorted(sigma_labels)[0])
    second_sigma = next(label for label in sorted(sigma_labels) if label != first_sigma)
    table = {
        "E": (1.0, 1.0, 1.0, 1.0),
        "C2x": (1.0, 1.0, -1.0, -1.0),
        "C2y": (1.0, 1.0, -1.0, -1.0),
        "C2z": (1.0, 1.0, -1.0, -1.0),
        "C2_perp": (1.0, 1.0, -1.0, -1.0),
        first_sigma: (1.0, -1.0, 1.0, -1.0),
        second_sigma: (1.0, -1.0, -1.0, 1.0),
    }
    arr = np.array([table[label] for label in labels], dtype=float)
    return [(name, arr[:, i]) for i, name in enumerate(("A1", "A2", "B1", "B2"))]


def _d2_irreps(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    return [
        ("A", _chars(labels, {"C2x": 1.0, "C2y": 1.0, "C2z": 1.0})),
        ("B1", _chars(labels, {"C2x": -1.0, "C2y": -1.0, "C2z": 1.0})),
        ("B2", _chars(labels, {"C2x": -1.0, "C2y": 1.0, "C2z": -1.0})),
        ("B3", _chars(labels, {"C2x": 1.0, "C2y": -1.0, "C2z": -1.0})),
    ]


def _c2h_irreps(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    c2 = {"C2x", "C2y", "C2z", "C2_perp"}
    return [
        ("Ag", np.array([1.0 if label == "E" or label in c2 or label == "i" or label.startswith("sigma") else 1.0 for label in labels])),
        ("Bg", np.array([1.0 if label in {"E", "i"} else -1.0 for label in labels])),
        ("Au", np.array([1.0 if label == "E" or label in c2 else -1.0 for label in labels])),
        ("Bu", np.array([1.0 if label == "E" or label.startswith("sigma") else -1.0 for label in labels])),
    ]


def _d2h_irreps(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    signs = {
        "Ag":  {"E": 1, "C2z": 1, "C2y": 1, "C2x": 1, "i": 1, "sigma_xy": 1, "sigma_xz": 1, "sigma_yz": 1},
        "B1g": {"E": 1, "C2z": 1, "C2y": -1, "C2x": -1, "i": 1, "sigma_xy": 1, "sigma_xz": -1, "sigma_yz": -1},
        "B2g": {"E": 1, "C2z": -1, "C2y": 1, "C2x": -1, "i": 1, "sigma_xy": -1, "sigma_xz": 1, "sigma_yz": -1},
        "B3g": {"E": 1, "C2z": -1, "C2y": -1, "C2x": 1, "i": 1, "sigma_xy": -1, "sigma_xz": -1, "sigma_yz": 1},
        "Au":  {"E": 1, "C2z": 1, "C2y": 1, "C2x": 1, "i": -1, "sigma_xy": -1, "sigma_xz": -1, "sigma_yz": -1},
        "B1u": {"E": 1, "C2z": 1, "C2y": -1, "C2x": -1, "i": -1, "sigma_xy": -1, "sigma_xz": 1, "sigma_yz": 1},
        "B2u": {"E": 1, "C2z": -1, "C2y": 1, "C2x": -1, "i": -1, "sigma_xy": 1, "sigma_xz": -1, "sigma_yz": 1},
        "B3u": {"E": 1, "C2z": -1, "C2y": -1, "C2x": 1, "i": -1, "sigma_xy": 1, "sigma_xz": 1, "sigma_yz": -1},
    }
    return [(name, _chars(labels, chars)) for name, chars in signs.items()]


def _cyclic_family_irreps(group: str, labels: list[str]) -> list[tuple[str, np.ndarray]]:
    match = re.match(r"C(\d+)$", group)
    if not match:
        return []
    n = int(match.group(1))
    out: list[tuple[str, np.ndarray]] = [("A", np.ones(len(labels), dtype=float))]
    if n % 2 == 0:
        out.append(("B", np.array([(-1.0) ** _rotation_power(label) for label in labels], dtype=float)))
    for k in range(1, (n + 1) // 2):
        if n % 2 == 0 and k == n // 2:
            continue
        chars = []
        for label in labels:
            chars.append(2.0 * np.cos(2.0 * np.pi * k * _rotation_power(label) / float(n)))
        out.append((f"E{k}", 2.0 * np.array(chars, dtype=float)))
    return out


def _rotation_power(label: str) -> int:
    if label == "E":
        return 0
    match = re.match(r"C(\d+)[xyz]\^(\d+)", label)
    if match:
        return int(match.group(2))
    if label.startswith("C2"):
        return 1
    return 0


def _oop_sign(orig, mapped):
    i, j, k, l = orig
    if mapped[1] != j:
        return 1.0
    outer = [i, k, l]
    mapped_outer = [mapped[0], mapped[2], mapped[3]]
    perm = [outer.index(x) if x in outer else -1 for x in mapped_outer]
    if -1 in perm:
        return 1.0
    inv = 0
    for a in range(3):
        for b in range(a + 1, 3):
            if perm[a] > perm[b]:
                inv += 1
    return -1.0 if (inv % 2) == 1 else 1.0


def primitive_permutation(prims, atom_map):
    n = len(prims)
    perm_idx = [-1] * n
    sign = [1.0] * n

    bond_map = {}
    angle_map = {}
    dihed_map = {}
    oop_map = {}
    linear_map = {}
    frag_map = {}

    for i, p in enumerate(prims):
        if p.kind == "bond":
            a, b = p.atoms
            bond_map.setdefault(tuple(sorted((a, b))), i)
        elif p.kind == "angle":
            a, j, b = p.atoms
            angle_map.setdefault((j, tuple(sorted((a, b)))), i)
        elif p.kind == "dihedral":
            dihed_map.setdefault(tuple(p.atoms), i)
        elif p.kind == "out_of_plane":
            oop_map.setdefault(tuple(p.atoms), i)
        elif p.kind == "linear_bend":
            a, j, b = p.atoms
            linear_map.setdefault((j, tuple(sorted((a, b))), p.mode), i)
        elif p.kind.startswith("frag_"):
            frag_map.setdefault((p.kind, tuple(sorted(p.atoms)), p.mode, tuple(sorted(p.ref))), i)

    for i, p in enumerate(prims):
        mapped_atoms = tuple(atom_map[a] for a in p.atoms)
        if p.kind == "bond":
            idx = bond_map.get(tuple(sorted(mapped_atoms)))
            if idx is not None:
                perm_idx[i] = idx
        elif p.kind == "angle":
            a, j, b = mapped_atoms
            idx = angle_map.get((j, tuple(sorted((a, b)))))
            if idx is not None:
                perm_idx[i] = idx
        elif p.kind == "dihedral":
            idx = dihed_map.get(mapped_atoms)
            if idx is None:
                idx = dihed_map.get(mapped_atoms[::-1])
                if idx is not None:
                    sign[i] = -1.0
            if idx is not None:
                perm_idx[i] = idx
        elif p.kind == "out_of_plane":
            idx = oop_map.get(mapped_atoms)
            if idx is None:
                for cand, ci in oop_map.items():
                    if cand[1] != mapped_atoms[1]:
                        continue
                    if set(cand) == set(mapped_atoms):
                        perm_idx[i] = ci
                        sign[i] = _oop_sign(p.atoms, cand)
                        break
            else:
                perm_idx[i] = idx
        elif p.kind == "linear_bend":
            a, j, b = mapped_atoms
            idx = linear_map.get((j, tuple(sorted((a, b))), p.mode))
            if idx is not None:
                perm_idx[i] = idx
        elif p.kind.startswith("frag_"):
            idx = frag_map.get((p.kind, tuple(sorted(mapped_atoms)), p.mode, tuple(sorted(p.ref))))
            if idx is not None:
                perm_idx[i] = idx

        if perm_idx[i] < 0:
            perm_idx[i] = i
            sign[i] = 1.0

    return perm_idx, sign


def symmetrize_u(
    U,
    prims,
    Z,
    coords,
    tol=1.0e-3,
    max_n=6,
    keep_a1_only=False,
    a1_tol=1e-6,
    per_type=True,
    label_symmetry=False,
    quasi_tol=None,
    tol_H=None,
    heavy_only_orient=False,
    center_idx=None,
    ignore_isotopes=False,
    op_filter=None,
    symbols_override=None,
    max_dev_strict=None,
    tol_rel=0.0,
    auto_max_n=False,
    inertia_tol=1e-3,
    max_radius=None,
    enforce_radial_filter=True,
    profile=False,
):
    symbols = list(symbols_override) if symbols_override is not None else [atomic_symbol(int(z)) for z in Z]
    weights = Z.copy()
    if heavy_only_orient:
        heavy = np.array([z for z in Z if z > 1], dtype=float)
        if len(heavy) >= 2:
            weights = np.array([z if z > 1 else 0.1 for z in Z], dtype=float)
    coords_oriented = orient_coords(coords, weights=weights)
    perf = {} if profile else None
    elements, classes, permutations = symmetry_elements_from_geometry(
        symbols,
        coords_oriented,
        tol=tol,
        max_n=max_n,
        tol_H=tol_H,
        ignore_isotopes=ignore_isotopes,
        op_filter=op_filter,
        tol_rel=tol_rel,
        auto_max_n=auto_max_n,
        inertia_tol=inertia_tol,
        max_radius=max_radius,
        enforce_radial_filter=enforce_radial_filter,
        perf=perf,
    )
    if max_dev_strict is not None and elements:
        filtered = [
            (elem, perm)
            for elem, perm in zip(elements, permutations)
            if elem[2] <= max_dev_strict
        ]
        if filtered:
            elements, permutations = zip(*filtered)
            elements = list(elements)
            permutations = list(permutations)
        else:
            elements, permutations = [], []

    linear = is_linear(coords_oriented, tol=tol)
    spherical_guess = _spherical_top_guess(
        symbols, coords_oriented, center_idx=center_idx, ignore_isotopes=ignore_isotopes
    )
    quasi_elements = None
    quasi_perms = None
    if quasi_tol is not None:
        quasi_elements, _, quasi_perms = symmetry_elements_from_geometry(
            symbols,
            coords_oriented,
            tol=quasi_tol,
            max_n=max_n,
            tol_H=tol_H,
            ignore_isotopes=ignore_isotopes,
            op_filter=op_filter,
            tol_rel=tol_rel,
            auto_max_n=auto_max_n,
            inertia_tol=inertia_tol,
            max_radius=max_radius,
            enforce_radial_filter=enforce_radial_filter,
            perf=perf,
        )
        if max_dev_strict is not None and quasi_elements:
            if quasi_perms is not None:
                filtered = [
                    (elem, perm)
                    for elem, perm in zip(quasi_elements, quasi_perms)
                    if elem[2] <= max_dev_strict
                ]
                if filtered:
                    quasi_elements, quasi_perms = zip(*filtered)
                    quasi_elements = list(quasi_elements)
                    quasi_perms = list(quasi_perms)
                else:
                    quasi_elements, quasi_perms = [], []
            else:
                quasi_elements = [e for e in quasi_elements if e[2] <= max_dev_strict]

    if not permutations:
        info = {
            "elements": [],
            "classes": classes,
            "point_group": "C1",
            "max_dev": 0.0,
            "mean_dev": 0.0,
            "generators": [],
            "nmax": 1,
            "axis": None,
        }
        return U, info

    perms = []
    for mapping in permutations:
        perm_idx, sign = primitive_permutation(prims, mapping)
        perms.append((perm_idx, sign))

    nprim, ncol = U.shape
    U_sym = np.zeros_like(U)
    if per_type:
        kinds = {}
        for i, p in enumerate(prims):
            kinds.setdefault(p.kind, []).append(i)
        for idxs in kinds.values():
            idxs_arr = np.array(idxs, dtype=int)
            loc = {gi: li for li, gi in enumerate(idxs)}
            for j in range(ncol):
                v = U[idxs_arr, j]
                acc = np.zeros(len(idxs_arr), dtype=float)
                for perm_idx, sign in perms:
                    perm_sub = [perm_idx[i] for i in idxs]
                    sign_sub = [sign[i] for i in idxs]
                    local_idx = [loc[p] for p in perm_sub]
                    vv = v[np.array(local_idx, dtype=int)] * np.array(sign_sub, dtype=float)
                    acc += vv
                U_sym[idxs_arr, j] = acc / len(perms)
    else:
        for j in range(ncol):
            v = U[:, j]
            acc = np.zeros(nprim, dtype=float)
            for perm_idx, sign in perms:
                vv = v[np.array(perm_idx, dtype=int)] * np.array(sign, dtype=float)
                acc += vv
            U_sym[:, j] = acc / len(perms)

    labels = None
    if label_symmetry:
        labels = []
        op_labels = [e[0] for e in elements]
        for j in range(U_sym.shape[1]):
            v = U_sym[:, j]
            vn = np.linalg.norm(v)
            if vn < 1e-12:
                labels.append({"label": "zero", "ops": {}})
                continue
            sig = {}
            all_pos = True
            for (perm_idx, sign), oplab in zip(perms, op_labels):
                vv = v[np.array(perm_idx, dtype=int)] * np.array(sign, dtype=float)
                c = float(np.dot(v, vv) / (vn * vn))
                if abs(c - 1.0) < a1_tol:
                    sig[oplab] = 1.0
                elif abs(c + 1.0) < a1_tol:
                    sig[oplab] = -1.0
                    all_pos = False
                else:
                    sig[oplab] = 0.0
                    all_pos = False
            labels.append({"label": "A1" if all_pos else "other", "ops": sig})

    if keep_a1_only:
        cols = []
        keep_labels = []
        for j in range(U_sym.shape[1]):
            if np.linalg.norm(U_sym[:, j]) > a1_tol:
                cols.append(U_sym[:, j])
                if labels is not None:
                    keep_labels.append(labels[j])
        if cols:
            U_sym = np.stack(cols, axis=1)
            if labels is not None:
                labels = keep_labels
        else:
            U_sym = np.zeros((nprim, 0), dtype=float)
            if labels is not None:
                labels = []

    max_dev = float(max((e[2] for e in elements), default=0.0))
    mean_dev = float(np.mean([e[2] for e in elements])) if elements else 0.0
    nmax, axis = _highest_cn_axis([e[0] for e in elements])
    info = {
        "elements": [e[0] for e in elements],
        "classes": classes,
        "point_group": spherical_guess if spherical_guess is not None else _group_label(elements, linear=linear),
        "max_dev": max_dev,
        "mean_dev": mean_dev,
        "generators": _pick_generators(elements),
        "nmax": nmax,
        "axis": axis,
    }
    if quasi_elements is not None:
        info["quasi_elements"] = [e[0] for e in quasi_elements]
        info["quasi_point_group"] = _group_label(quasi_elements, linear=linear)
        info["quasi_max_dev"] = float(max((e[2] for e in quasi_elements), default=0.0))
    if labels is not None:
        info["labels"] = labels
    if perf is not None:
        info["perf"] = perf
    return U_sym, info
