from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from merlino_fit.survibfit.modify_geom import read_xyz
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
    op_data = _operation_data(atoms, coords, prims)
    sym_gics = _symmetry_adapted_gics(gics, prims, u_matrix, op_data)
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


def _operation_data(atoms: list[str], coords: np.ndarray, prims: list[Primitive]):
    z_numbers = np.array([atomic_number(atom) for atom in atoms], dtype=int)
    symbols = [atomic_symbol(int(z)) for z in z_numbers]
    oriented = orient_coords(coords, weights=z_numbers)
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
        unique.append((element[0], mapped, primitive_permutation(prims, mapped)))
    return unique or [("E", tuple(range(len(atoms))), primitive_permutation(prims, tuple(range(len(atoms)))))]


def _symmetry_adapted_gics(gics, prims, u_matrix, op_data):
    irreps = _irrep_characters([item[0] for item in op_data])
    if not irreps:
        return [(gic.name, "A", u_matrix[:, idx]) for idx, gic in enumerate(gics)]
    adapted = []
    used_names: dict[str, int] = {}
    for kind in _kind_order(prims):
        idxs = np.array([i for i, prim in enumerate(prims) if prim.kind == kind], dtype=int)
        block = u_matrix[idxs, :]
        for irrep, chars in irreps:
            basis: list[np.ndarray] = []
            for col in range(block.shape[1]):
                projected = np.zeros(block.shape[0], dtype=float)
                source = u_matrix[:, col]
                for op_index, (_label, _mapping, (perm_idx, sign)) in enumerate(op_data):
                    transformed = np.zeros_like(source)
                    for src, dst in enumerate(perm_idx):
                        transformed[dst] += sign[src] * source[src]
                    projected += chars[op_index] * transformed[idxs]
                projected /= float(len(op_data))
                if np.linalg.norm(projected) < 1.0e-8:
                    continue
                projected = _orthogonal_residual(projected, basis)
                norm = np.linalg.norm(projected)
                if norm < 1.0e-7:
                    continue
                basis.append(projected / norm)
                full = np.zeros(len(prims), dtype=float)
                full[idxs] = basis[-1]
                adapted.append((_next_name(irrep, kind, used_names), irrep, full))
    return adapted


def _kind_order(prims: list[Primitive]) -> list[str]:
    order = ["bond", "angle", "linear_bend", "dihedral", "out_of_plane"]
    present = {prim.kind for prim in prims}
    return [kind for kind in order if kind in present]


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
    for name, _irrep, column in sym_gics:
        out.append(_format_gic_line(name, column, prims))
    target.write_text("\n".join(out) + "\n", encoding="utf-8")


def _format_gic_line(name: str, column: np.ndarray, prims: list[Primitive]) -> str:
    parts = []
    for coeff, primitive in zip(column, prims):
        if abs(coeff) < 1.0e-8:
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
