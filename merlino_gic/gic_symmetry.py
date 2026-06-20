from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from merlino_fit.survibfit.modify_geom import read_xyz
from merlino_fit.survibfit.pipeline import b_matrix_analytic
from merlino_fit.survibfit.primitives import Primitive
from merlino_fit.survibfit.symmetry_detector import orient_coords, symmetry_elements_from_geometry
from merlino_fit.survibfit.symmetry_global import primitive_permutation
from topology.elements import atomic_number, atomic_symbol


@dataclass(frozen=True)
class GICLine:
    name: str
    terms: tuple[tuple[float, Primitive], ...]


def write_gic_symmetry_files(workdir: Path) -> None:
    run_dir = Path(workdir)
    gauin = run_dir / "gauin"
    xyzin = run_dir / "xyzin"
    if not gauin.exists() or not xyzin.exists():
        return
    atoms, coords, _comment = read_xyz(xyzin)
    gics = _parse_gauin_gics(gauin)
    if not gics:
        return
    prims, u_matrix = _primitive_basis(gics)
    oriented = _oriented_coords(atoms, coords)
    op_data = _operation_data(atoms, oriented, prims, already_oriented=True)
    sym_gics = _symmetry_adapted_gics(gics, prims, u_matrix, op_data, oriented)
    _write_gicsym(run_dir / "gicsym", sym_gics)
    _write_symmetrized_gauin(gauin, run_dir / "gauin.symm", sym_gics, prims)


def _parse_gauin_gics(gauin: Path) -> list[GICLine]:
    out: list[GICLine] = []
    for raw in gauin.read_text(encoding="utf-8", errors="replace").splitlines():
        parsed = _parse_gic_line(raw)
        if parsed is not None:
            out.append(parsed)
    return out


def _parse_gic_line(line: str) -> GICLine | None:
    stripped = line.strip()
    if not stripped or "=" not in stripped:
        return None
    name, rhs = stripped.split("=", 1)
    name = name.replace("(Inactive)", "").strip()
    if name.startswith(("QPck", "PhiP")):
        return None
    number = r"[+-]?\s*(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?"
    terms: list[tuple[float, Primitive]] = []
    for match in re.finditer(rf"({number})\s*\*\s*([RADLU])\(([^)]*)\)", rhs):
        coeff = float(match.group(1).replace(" ", "").replace("D", "E").replace("d", "e"))
        terms.append((coeff, _primitive(match.group(2), match.group(3))))
    if not terms:
        simple = re.search(r"\b([RADLU])\(([^)]*)\)", rhs)
        if simple:
            terms.append((1.0, _primitive(simple.group(1), simple.group(2))))
    if not terms:
        return None
    return GICLine(name=name, terms=tuple(terms))


def _primitive(kind: str, atoms_text: str) -> Primitive:
    values = tuple(int(item.strip()) for item in atoms_text.split(",") if item.strip())
    atoms = tuple(value - 1 for value in values)
    if kind == "R" and len(atoms) == 2:
        return Primitive("bond", atoms)
    if kind == "A" and len(atoms) == 3:
        return Primitive("angle", atoms)
    if kind == "D" and len(atoms) == 4:
        return Primitive("dihedral", atoms)
    if kind == "U" and len(atoms) == 4:
        return Primitive("out_of_plane", atoms)
    if kind == "L" and len(values) == 5:
        mode = values[4] if values[4] in {-1, -2} else -1
        return Primitive("linear_bend", atoms[:3], mode=mode)
    raise ValueError(f"Unsupported GIC primitive {kind}({atoms_text})")


def _primitive_basis(gics: list[GICLine]) -> tuple[list[Primitive], np.ndarray]:
    prims: list[Primitive] = []
    index: dict[Primitive, int] = {}
    columns: list[np.ndarray] = []
    for gic in gics:
        col = np.zeros(len(prims), dtype=float)
        for coeff, prim in gic.terms:
            if prim not in index:
                index[prim] = len(prims)
                prims.append(prim)
                col = np.pad(col, (0, 1))
                for i, existing in enumerate(columns):
                    columns[i] = np.pad(existing, (0, 1))
            col[index[prim]] += coeff
        columns.append(col)
    return prims, np.column_stack(columns)


def _oriented_coords(atoms: list[str], coords: np.ndarray) -> np.ndarray:
    z_numbers = np.array([atomic_number(atom) for atom in atoms], dtype=int)
    return orient_coords(coords, weights=z_numbers)


def _operation_data(atoms: list[str], coords: np.ndarray, prims: list[Primitive], already_oriented: bool = False):
    z_numbers = np.array([atomic_number(atom) for atom in atoms], dtype=int)
    symbols = [atomic_symbol(int(z)) for z in z_numbers]
    oriented = coords if already_oriented else orient_coords(coords, weights=z_numbers)
    elements, _classes, permutations = symmetry_elements_from_geometry(
        symbols,
        oriented,
        tol=1.0e-2,
        max_n=6,
        tol_H=1.0e-2,
        ignore_isotopes=True,
        auto_max_n=True,
        inertia_tol=1.0e-3,
    )
    unique = []
    seen = set()
    for element, permutation in zip(elements, permutations):
        mapped = tuple(int(item) for item in permutation)
        if mapped in seen:
            continue
        seen.add(mapped)
        unique.append((element[0], element[1], mapped, primitive_permutation(prims, mapped)))
    identity = tuple(range(len(atoms)))
    return unique or [("E", np.eye(3), identity, primitive_permutation(prims, identity))]


def _symmetry_adapted_gics(gics, prims, u_matrix, op_data, coords: np.ndarray):
    irreps = _irrep_characters([item[0] for item in op_data])
    if not irreps:
        return [(gic.name, "A", u_matrix[:, idx]) for idx, gic in enumerate(gics)]
    targets = _vibrational_irrep_counts(op_data, irreps, len(coords))
    b_primitive = b_matrix_analytic(prims, coords)
    source_rows = u_matrix.T @ b_primitive
    vib_projector = _vibrational_projector(coords)
    cart_ops = [_cartesian_operation(rotation, mapping, len(coords)) for _label, rotation, mapping, _prim_op in op_data]
    adapted = []
    used_names: dict[str, int] = {}
    selected_rows: dict[str, list[np.ndarray]] = {irrep: [] for irrep, _chars in irreps}
    selected_global: list[np.ndarray] = []
    for irrep, chars in irreps:
        target = targets.get(irrep, 0)
        if target <= 0:
            continue
        for col, source_row in enumerate(source_rows):
            projected_row_raw = _project_cartesian_row(source_row, chars, cart_ops)
            projected_row = projected_row_raw @ vib_projector
            if np.linalg.norm(projected_row) < 1.0e-8:
                continue
            residual = _orthogonal_residual(projected_row, selected_rows[irrep])
            row_norm = np.linalg.norm(residual)
            if row_norm < 1.0e-7:
                continue
            global_residual = _orthogonal_residual(projected_row, selected_global)
            if np.linalg.norm(global_residual) < 1.0e-7:
                continue
            coeff = _project_column_to_irrep(u_matrix[:, col], chars, op_data)
            coeff_norm = np.linalg.norm(coeff)
            if coeff_norm < 1.0e-8:
                coeff = _cartesian_row_to_coeff(projected_row_raw, b_primitive, u_matrix[:, col], prims)
                coeff_norm = np.linalg.norm(coeff)
                if coeff_norm < 1.0e-8:
                    continue
            coeff /= coeff_norm
            selected_rows[irrep].append(residual / row_norm)
            selected_global.append(global_residual / np.linalg.norm(global_residual))
            kind = _dominant_kind(u_matrix[:, col], prims)
            adapted.append((_next_name(irrep, kind, used_names), irrep, coeff))
            if len(selected_rows[irrep]) == target:
                break
        if len(selected_rows[irrep]) != target:
            raise RuntimeError(
                f"GIC symmetry reduction generated {len(selected_rows[irrep])} {irrep} coordinates; expected {target}"
            )
    counts = {irrep: len(rows) for irrep, rows in selected_rows.items()}
    if counts != targets:
        raise RuntimeError(f"GIC symmetry reduction count mismatch: {counts}; expected {targets}")
    return adapted


def _cartesian_operation(rotation: np.ndarray, mapping: tuple[int, ...], natoms: int) -> np.ndarray:
    matrix = np.zeros((3 * natoms, 3 * natoms), dtype=float)
    # The detector returns i -> j such that x_i matches R x_j; for row
    # gradients this block form applies the same operation in the oriented
    # Cartesian frame.
    for i, j in enumerate(mapping):
        matrix[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = rotation
    return matrix


def _project_cartesian_row(row: np.ndarray, chars: np.ndarray, cart_ops: list[np.ndarray]) -> np.ndarray:
    projected = np.zeros_like(row)
    for op_index, op_matrix in enumerate(cart_ops):
        projected += chars[op_index] * (row @ op_matrix)
    return projected / float(len(cart_ops))


def _vibrational_projector(coords: np.ndarray) -> np.ndarray:
    natoms = len(coords)
    basis = []
    for axis in range(3):
        vec = np.zeros(3 * natoms, dtype=float)
        vec[axis::3] = 1.0
        basis.append(vec)
    for axis in np.eye(3):
        vec = np.array([component for coord in coords for component in np.cross(axis, coord)], dtype=float)
        basis.append(vec)
    ortho: list[np.ndarray] = []
    for vec in basis:
        residual = _orthogonal_residual(vec, ortho)
        norm = np.linalg.norm(residual)
        if norm > 1.0e-10:
            ortho.append(residual / norm)
    if not ortho:
        return np.eye(3 * natoms, dtype=float)
    q_matrix = np.vstack(ortho).T
    return np.eye(3 * natoms, dtype=float) - q_matrix @ q_matrix.T


def _cartesian_row_to_coeff(
    row: np.ndarray, b_primitive: np.ndarray, source_coeff: np.ndarray, prims: list[Primitive]
) -> np.ndarray:
    kind = _dominant_kind(source_coeff, prims)
    idxs = [idx for idx, prim in enumerate(prims) if prim.kind == kind]
    coeff = _least_squares_coeff(row, b_primitive, idxs)
    if np.linalg.norm(coeff @ b_primitive - row) > 1.0e-5 * max(1.0, np.linalg.norm(row)):
        coeff = _least_squares_coeff(row, b_primitive, list(range(len(prims))))
    return coeff


def _least_squares_coeff(row: np.ndarray, b_primitive: np.ndarray, idxs: list[int]) -> np.ndarray:
    coeff = np.zeros(b_primitive.shape[0], dtype=float)
    if not idxs:
        return coeff
    sub_b = b_primitive[np.array(idxs, dtype=int), :]
    values, *_ = np.linalg.lstsq(sub_b.T, row.T, rcond=1.0e-10)
    coeff[np.array(idxs, dtype=int)] = values
    return coeff


def _project_column_to_irrep(source: np.ndarray, chars: np.ndarray, op_data) -> np.ndarray:
    projected = np.zeros_like(source)
    for op_index, (_label, _rotation, _mapping, (perm_idx, sign)) in enumerate(op_data):
        transformed = np.zeros_like(source)
        for src, dst in enumerate(perm_idx):
            transformed[dst] += sign[src] * source[src]
        projected += chars[op_index] * transformed
    return projected / float(len(op_data))


def _dominant_kind(column: np.ndarray, prims: list[Primitive]) -> str:
    weights: dict[str, float] = {}
    for coeff, prim in zip(column, prims):
        weights[prim.kind] = weights.get(prim.kind, 0.0) + float(coeff * coeff)
    if not weights:
        return "gic"
    return max(weights.items(), key=lambda item: item[1])[0]


def _orthogonal_residual(vector: np.ndarray, basis: list[np.ndarray]) -> np.ndarray:
    residual = vector.astype(float, copy=True)
    for item in basis:
        residual -= np.dot(item, residual) * item
    return residual


def _irrep_characters(labels: list[str]) -> list[tuple[str, np.ndarray]]:
    if len(labels) == 1:
        return [("A", np.ones(1))]
    if len(labels) == 2:
        return [("A", np.array([1.0, 1.0])), ("B", np.array([1.0, -1.0]))]
    if len(labels) == 4 and any(label.startswith("C2") for label in labels) and sum(label.startswith("sigma") for label in labels) == 2:
        chars = []
        for label in labels:
            if label == "E":
                chars.append((1.0, 1.0, 1.0, 1.0))
            elif label.startswith("C2"):
                chars.append((1.0, 1.0, -1.0, -1.0))
            elif label == "sigma_xz":
                chars.append((1.0, -1.0, 1.0, -1.0))
            else:
                chars.append((1.0, -1.0, -1.0, 1.0))
        arr = np.array(chars, dtype=float)
        return [(name, arr[:, i]) for i, name in enumerate(("A1", "A2", "B1", "B2"))]
    return []


def _vibrational_irrep_counts(op_data, irreps: list[tuple[str, np.ndarray]], natoms: int) -> dict[str, int]:
    gamma_3n = []
    gamma_trans = []
    gamma_rot = []
    for _label, rotation, mapping, _primitive_op in op_data:
        fixed = sum(1 for i, j in enumerate(mapping) if i == j)
        trace = float(np.trace(rotation))
        gamma_3n.append(fixed * trace)
        gamma_trans.append(trace)
        gamma_rot.append(float(np.linalg.det(rotation) * trace))
    gamma_vib = np.array(gamma_3n) - np.array(gamma_trans) - np.array(gamma_rot)
    counts = {}
    group_order = float(len(op_data))
    for irrep, chars in irreps:
        value = int(round(float(np.dot(gamma_vib, chars)) / group_order))
        counts[irrep] = max(value, 0)
    if sum(counts.values()) != 3 * natoms - 6:
        raise RuntimeError(
            f"Vibrational irrep count mismatch: {counts} sums to {sum(counts.values())}, expected {3 * natoms - 6}"
        )
    return counts


def _next_name(irrep: str, kind: str, used: dict[str, int]) -> str:
    prefix = {
        "bond": "Str",
        "angle": "Ang",
        "linear_bend": "Lin",
        "dihedral": "Tor",
        "out_of_plane": "Oop",
    }.get(kind, "Gic")
    key = f"{irrep}{prefix}"
    used[key] = used.get(key, 0) + 1
    return f"{irrep}{prefix}{used[key]:04d}"


def _write_gicsym(path: Path, sym_gics) -> None:
    lines = ["name,irrep"]
    for name, irrep, _column in sym_gics:
        lines.append(f"{name},{irrep}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_symmetrized_gauin(source: Path, target: Path, sym_gics, prims: list[Primitive]) -> None:
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    first_gic = next((i for i, line in enumerate(lines) if _parse_gic_line(line) is not None), len(lines))
    prefix = lines[:first_gic]
    out = list(prefix)
    a1 = [item for item in sym_gics if item[1] in {"A1", "A", "Ag", "A'"}]
    other = [item for item in sym_gics if item not in a1]
    for name, _irrep, column in a1:
        out.append(_format_gic_line(name, column, prims))
    if other:
        out.append("")
    for name, _irrep, column in other:
        out.append(_format_gic_line(name, column, prims))
    target.write_text("\n".join(out) + "\n", encoding="utf-8")


def _format_gic_line(name: str, column: np.ndarray, prims: list[Primitive]) -> str:
    parts = []
    for coeff, primitive in zip(column, prims):
        if abs(coeff) < 1.0e-5:
            continue
        parts.append((coeff, _primitive_expression(primitive)))
    expr = _join_terms(parts)
    return f" {name}=[ {expr}]"


def _primitive_expression(primitive: Primitive) -> str:
    atoms = tuple(atom + 1 for atom in primitive.atoms)
    if primitive.kind == "bond":
        return f"R({atoms[0]:3d},{atoms[1]:3d})"
    if primitive.kind == "angle":
        return f"A({atoms[0]:3d},{atoms[1]:3d},{atoms[2]:3d})"
    if primitive.kind == "dihedral":
        return f"D({atoms[0]:3d},{atoms[1]:3d},{atoms[2]:3d},{atoms[3]:3d})"
    if primitive.kind == "out_of_plane":
        return f"U({atoms[0]:3d},{atoms[1]:3d},{atoms[2]:3d},{atoms[3]:3d})"
    if primitive.kind == "linear_bend":
        return f"L({atoms[0]:3d},{atoms[1]:3d},{atoms[2]:3d},  0,{primitive.mode:3d})"
    raise ValueError(f"Unsupported primitive kind: {primitive.kind}")


def _join_terms(parts: list[tuple[float, str]]) -> str:
    text = ""
    for idx, (coeff, expr) in enumerate(parts):
        sign = "-" if coeff < 0.0 else "+"
        body = f"{abs(coeff):.8f}*{expr}"
        if idx == 0:
            text += f"-{body}" if sign == "-" else body
        else:
            text += sign + body
    return text
