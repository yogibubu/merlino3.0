from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from merlino_fit.survibfit.geometry import angle
from merlino_fit.survibfit.pipeline import b_matrix_analytic
from merlino_fit.survibfit.primitives import Primitive
from merlino_fit.topology.pipeline import build_topology_objects
from topology.elements import atomic_number

from .model import (
    GICDefinition,
    _definition_coordinate_kind_counts,
    _gicforge_cartesian_from_gauin,
    _primitive_signature,
    define_gics_from_cartesian,
)


LINEAR_THRESHOLD_RAD = np.deg2rad(170.0)


@dataclass(frozen=True)
class GICForgePythonCoordinate:
    name: str
    block: str
    terms: tuple[tuple[float, Primitive], ...]
    type_index: int = 0

    @property
    def dominant_kind(self) -> str:
        return self.terms[0][1].kind if self.terms else "unknown"


@dataclass(frozen=True)
class GICForgePythonModel:
    atom_symbols: tuple[str, ...]
    atomic_numbers: tuple[int, ...]
    coordinates_angstrom: tuple[tuple[float, float, float], ...]
    primitive_candidates: tuple[GICForgePythonCoordinate, ...]
    coordinates: tuple[GICForgePythonCoordinate, ...]
    target_rank: int
    primitive_fallback: bool

    def to_definition(self, *, workdir: Path | None = None) -> GICDefinition:
        primitive_basis = _primitive_basis(self.coordinates)
        row_index = {primitive: index for index, primitive in enumerate(primitive_basis)}
        u_matrix = np.zeros((len(primitive_basis), len(self.coordinates)), dtype=float)
        labels: list[str] = []
        names: list[str] = []
        for column, coordinate in enumerate(self.coordinates):
            names.append(coordinate.name)
            for coefficient, primitive in coordinate.terms:
                u_matrix[row_index[primitive], column] += float(coefficient)
            labels.append(
                f"GIC{column + 1:03d} GICForgePython {coordinate.name} "
                f"irrep=UNK {_format_terms(coordinate.terms)}"
            )
        definition = GICDefinition(
            atom_symbols=self.atom_symbols,
            atomic_numbers=self.atomic_numbers,
            reference_coordinates_angstrom=self.coordinates_angstrom,
            primitives=primitive_basis,
            u_matrix=u_matrix,
            labels=tuple(labels),
            names=tuple(names),
            irreps=tuple("UNK" for _ in names),
            point_group="UNKNOWN",
            symmetrized=False,
            symmetry_source="none",
            gaussian_input="\n".join(_format_readgic(name, coord.terms) for name, coord in zip(names, self.coordinates))
            + "\n",
            source="gicforge-python",
            generation_workdir=str(workdir) if workdir is not None else None,
            provenance={
                "backend": "gicforge-python",
                "target_vibrational_rank": str(self.target_rank),
                "primitive_fallback": str(self.primitive_fallback).lower(),
            },
        )
        return definition


def build_gicforge_python_model(
    atom_symbols: Iterable[str],
    coordinates_angstrom: np.ndarray,
    *,
    impdih: bool = True,
    linear_threshold: float = LINEAR_THRESHOLD_RAD,
    primitive_fallback: bool = True,
) -> GICForgePythonModel:
    atoms = tuple(str(atom).strip() for atom in atom_symbols)
    coords = np.asarray(coordinates_angstrom, dtype=float)
    if coords.shape != (len(atoms), 3):
        raise ValueError(f"Expected coordinate shape ({len(atoms)}, 3), got {coords.shape}")
    atomic_numbers = tuple(atomic_number(atom) for atom in atoms)
    _cg, graph, ringset, _synthons, _aromaticity = build_topology_objects(coords, np.asarray(atomic_numbers))
    primitive_blocks = _fortran_like_primitive_blocks(
        graph,
        coords,
        atomic_numbers=atomic_numbers,
        ringset=ringset,
        impdih=impdih,
        linear_threshold=linear_threshold,
    )
    primitive_candidates = tuple(coord for block in primitive_blocks for coord in block)
    target = _target_rank(coords, graph)
    candidates = primitive_candidates
    if not primitive_fallback and len(candidates) < target:
        raise ValueError(f"GICForge Python candidates below vibrational rank ({len(candidates)} < {target})")
    if len(candidates) < target:
        raise ValueError(f"Primitive candidates below vibrational rank ({len(candidates)} < {target})")
    coordinates = _prune_type_local(candidates, coords, target_rank=target)
    return GICForgePythonModel(
        atom_symbols=atoms,
        atomic_numbers=atomic_numbers,
        coordinates_angstrom=tuple(tuple(float(value) for value in row) for row in coords),
        primitive_candidates=primitive_candidates,
        coordinates=coordinates,
        target_rank=target,
        primitive_fallback=True,
    )


def compare_gicforge_python_to_fortran(
    atom_symbols: Iterable[str],
    coordinates_angstrom: np.ndarray,
    *,
    workdir: Path,
    executable: Path | None = None,
    impdih: bool = True,
) -> dict[str, object]:
    workdir = Path(workdir)
    fortran_dir = workdir / "fortran"
    python_model = build_gicforge_python_model(atom_symbols, coordinates_angstrom, impdih=impdih)
    fortran_definition = define_gics_from_cartesian(
        tuple(atom_symbols),
        np.asarray(coordinates_angstrom, dtype=float),
        workdir=fortran_dir,
        executable=executable,
        symmetrize=False,
    )
    raw_coords = _gicforge_cartesian_from_gauin(fortran_dir / "gauin", len(fortran_definition.atom_symbols))
    python_on_fortran_frame = build_gicforge_python_model(
        atom_symbols,
        raw_coords,
        impdih=impdih,
    ).to_definition(workdir=workdir / "python")
    fortran_signatures = tuple(_primitive_signature(primitive) for primitive in fortran_definition.primitives)
    python_signatures = tuple(_primitive_signature(primitive) for primitive in python_on_fortran_frame.primitives)
    same_ordered_primitives = fortran_signatures == python_signatures
    b_max_abs_diff = None
    if same_ordered_primitives and fortran_definition.u_matrix.shape == python_on_fortran_frame.u_matrix.shape:
        fortran_b = fortran_definition.u_matrix.T @ b_matrix_analytic(fortran_definition.primitives, raw_coords)
        python_b = python_on_fortran_frame.u_matrix.T @ b_matrix_analytic(
            python_on_fortran_frame.primitives,
            raw_coords,
        )
        b_max_abs_diff = float(np.max(np.abs(python_b - fortran_b))) if python_b.size else 0.0
    return {
        "passed": (
            len(python_on_fortran_frame.names) == len(fortran_definition.names)
            and _definition_coordinate_kind_counts(python_on_fortran_frame)
            == _definition_coordinate_kind_counts(fortran_definition)
            and same_ordered_primitives
            and (b_max_abs_diff is None or b_max_abs_diff <= 1.0e-7)
        ),
        "target_rank": python_model.target_rank,
        "python_gic_count": len(python_on_fortran_frame.names),
        "fortran_gic_count": len(fortran_definition.names),
        "python_kind_counts": _definition_coordinate_kind_counts(python_on_fortran_frame),
        "fortran_kind_counts": _definition_coordinate_kind_counts(fortran_definition),
        "python_primitive_count": len(python_on_fortran_frame.primitives),
        "fortran_primitive_count": len(fortran_definition.primitives),
        "same_ordered_primitives": same_ordered_primitives,
        "b_max_abs_diff": b_max_abs_diff,
        "python_names": list(python_on_fortran_frame.names),
        "fortran_names": list(fortran_definition.names),
        "python_workdir": str(workdir / "python"),
        "fortran_workdir": str(fortran_dir),
    }


def _fortran_like_primitive_blocks(
    graph,
    coords: np.ndarray,
    *,
    atomic_numbers: tuple[int, ...],
    ringset,
    impdih: bool,
    linear_threshold: float,
):
    bonds: list[GICForgePythonCoordinate] = []
    bends: list[GICForgePythonCoordinate] = []
    linears: list[GICForgePythonCoordinate] = []
    torsions: list[GICForgePythonCoordinate] = []
    oops: list[GICForgePythonCoordinate] = []
    neighbors = [sorted(graph.adjacency[index]) for index in range(graph.natoms)]
    atom_ring = _atom_ring_map(ringset, graph.natoms)

    for center in range(graph.natoms):
        neigh = neighbors[center]
        if len(neigh) == 3:
            bends.extend(
                _c2v3_angle_coordinates(
                    center,
                    neigh,
                    atomic_numbers=atomic_numbers,
                    neighbors=neighbors,
                    atom_ring=atom_ring,
                    start=len(bends) + 1,
                )
            )
        elif len(neigh) > 1:
            if len(neigh) > 3:
                # TODO: port FourAt/HighCoordAt exactly. For now keep primitive
                # angles so the contract reports the remaining gap explicitly.
                pass
            for ib, first_angle in enumerate(neigh[:-1]):
                for second_angle in neigh[ib + 1 :]:
                    value = angle(first_angle, center, second_angle, coords)
                    left, right = sorted((first_angle, second_angle))
                    if value < linear_threshold:
                        bends.append(
                            _primitive_coordinate("Bend", len(bends) + 1, Primitive("angle", (left, center, right)))
                        )
                    else:
                        linears.append(
                            _primitive_coordinate(
                                "LAng",
                                len(linears) + 1,
                                Primitive("linear_bend", (left, center, right), mode=-1),
                            )
                        )
                        linears.append(
                            _primitive_coordinate(
                                "LAng",
                                len(linears) + 1,
                                Primitive("linear_bend", (left, center, right), mode=-2),
                            )
                        )
        for ib, first in enumerate(neigh):
            if first < center:
                continue
            bonds.append(_primitive_coordinate("Stre", len(bonds) + 1, Primitive("bond", (center, first))))

    for bond in bonds:
        _coef, primitive = bond.terms[0]
        center, right = primitive.atoms
        if len(neighbors[center]) == 1 or len(neighbors[right]) == 1:
            continue
        for left in neighbors[center]:
            if left == right:
                continue
            if angle(left, center, right, coords) > linear_threshold:
                continue
            for far in neighbors[right]:
                if far == center:
                    continue
                if angle(center, right, far, coords) > linear_threshold:
                    continue
                if far == left:
                    continue
                torsions.append(
                    _primitive_coordinate("Tors", len(torsions) + 1, Primitive("dihedral", (left, center, right, far)))
                )

    oop_prefix = "ImpD" if impdih else "OuPl"
    for center in range(graph.natoms):
        neigh = neighbors[center]
        if len(neigh) != 3:
            continue
        first, second, third = neigh
        if impdih:
            primitive = Primitive("dihedral", (first, center, third, second))
        else:
            primitive = Primitive("out_of_plane", (center, first, second, third))
        oops.append(_primitive_coordinate(oop_prefix, len(oops) + 1, primitive))

    return bonds, bends, linears, torsions, oops


def _c2v3_angle_coordinates(
    center: int,
    neigh: list[int],
    *,
    atomic_numbers: tuple[int, ...],
    neighbors: list[list[int]],
    atom_ring: list[int],
    start: int,
) -> list[GICForgePythonCoordinate]:
    first, second, third = neigh
    eq12 = atomic_numbers[first] == atomic_numbers[second]
    eq13 = atomic_numbers[first] == atomic_numbers[third]
    eq23 = atomic_numbers[second] == atomic_numbers[third]
    if not eq12 and not eq13 and not eq23:
        different = first
        if atomic_numbers[second] == 1:
            different = second
        elif atomic_numbers[third] == 1:
            different = third
        elif len(neighbors[second]) == 1:
            different = second
        elif len(neighbors[third]) == 1:
            different = third
    elif eq12 and eq13:
        different = first
    elif eq12 and not eq13:
        different = third
    elif eq13 and not eq12:
        different = second
    else:
        different = first

    if atom_ring[center] != 0:
        if atom_ring[first] != 0 and atom_ring[second] != 0:
            different = third
        if atom_ring[first] != 0 and atom_ring[third] != 0:
            different = second
        if atom_ring[second] != 0 and atom_ring[third] != 0:
            different = first

    if different == first:
        jat, kat, lat = first, second, third
    elif different == second:
        jat, kat, lat = second, first, third
    else:
        jat, kat, lat = third, second, first

    if atom_ring[center] != 0 and atom_ring[jat] != 0:
        return []

    coords: list[GICForgePythonCoordinate] = []
    if atom_ring[center] == 0:
        den = np.sqrt(6.0)
        coords.append(
            GICForgePythonCoordinate(
                name=f"SymD{start:04d}",
                block="SymD",
                type_index=1,
                terms=(
                    (2.0 / den, Primitive("angle", tuple(sorted((kat, lat))[:1]) + (center,) + tuple(sorted((kat, lat))[1:]))),
                    (-1.0 / den, Primitive("angle", tuple(sorted((jat, kat))[:1]) + (center,) + tuple(sorted((jat, kat))[1:]))),
                    (-1.0 / den, Primitive("angle", tuple(sorted((jat, lat))[:1]) + (center,) + tuple(sorted((jat, lat))[1:]))),
                ),
            )
        )
        start += 1
    den = np.sqrt(2.0)
    coords.append(
        GICForgePythonCoordinate(
            name=f"Rock{start:04d}",
            block="Rock",
            type_index=2,
            terms=(
                (1.0 / den, Primitive("angle", tuple(sorted((jat, kat))[:1]) + (center,) + tuple(sorted((jat, kat))[1:]))),
                (-1.0 / den, Primitive("angle", tuple(sorted((jat, lat))[:1]) + (center,) + tuple(sorted((jat, lat))[1:]))),
            ),
        )
    )
    return coords


def _atom_ring_map(ringset, natoms: int) -> list[int]:
    atom_ring = [0 for _ in range(natoms)]
    if ringset is None:
        return atom_ring
    for index, ring in enumerate(ringset, start=1):
        for atom in ring.atoms:
            atom_ring[int(atom)] = index
    return atom_ring


def _primitive_coordinate(prefix: str, index: int, primitive: Primitive) -> GICForgePythonCoordinate:
    return GICForgePythonCoordinate(
        name=f"{prefix}{index:04d}",
        block=prefix,
        terms=((1.0, primitive),),
    )


def _prune_type_local(
    coordinates: tuple[GICForgePythonCoordinate, ...],
    coords: np.ndarray,
    *,
    target_rank: int,
) -> tuple[GICForgePythonCoordinate, ...]:
    by_kind = {
        "bond": [coord for coord in coordinates if coord.dominant_kind == "bond"],
        "angle": [coord for coord in coordinates if coord.dominant_kind == "angle"],
        "linear_bend": [coord for coord in coordinates if coord.dominant_kind == "linear_bend"],
        "dihedral": [coord for coord in coordinates if coord.dominant_kind == "dihedral"],
        "out_of_plane": [coord for coord in coordinates if coord.dominant_kind == "out_of_plane"],
    }
    ordered = by_kind["bond"] + by_kind["angle"] + by_kind["linear_bend"] + by_kind["dihedral"] + by_kind["out_of_plane"]
    if len(ordered) <= target_rank:
        return tuple(ordered)
    primitive_basis = _primitive_basis(ordered)
    row_index = {primitive: index for index, primitive in enumerate(primitive_basis)}
    primitive_b = b_matrix_analytic(primitive_basis, coords)
    b_rows = []
    for coordinate in ordered:
        row = np.zeros(primitive_b.shape[1], dtype=float)
        for coefficient, primitive in coordinate.terms:
            row += coefficient * primitive_b[row_index[primitive]]
        b_rows.append(row)

    basis: list[np.ndarray] = []
    keep: list[GICForgePythonCoordinate] = []
    for index, coordinate in enumerate(ordered):
        if coordinate.dominant_kind == "bond":
            _seed_basis_row(b_rows[index], basis)
            keep.append(coordinate)
    for kind in ("angle", "linear_bend", "dihedral", "out_of_plane"):
        for index, coordinate in enumerate(ordered):
            if coordinate.dominant_kind != kind:
                continue
            if len(basis) >= target_rank:
                continue
            if _seed_basis_row(b_rows[index], basis):
                keep.append(coordinate)
    return tuple(keep)


def _seed_basis_row(row: np.ndarray, basis: list[np.ndarray]) -> bool:
    t_abs = 1.0e-10
    t_rel = 1.0e-8
    candidate = np.asarray(row, dtype=float).copy()
    norm0 = float(np.linalg.norm(candidate))
    if norm0 <= t_abs:
        return False
    for existing in basis:
        candidate -= float(np.dot(candidate, existing)) * existing
    norm = float(np.linalg.norm(candidate))
    if norm > t_abs and norm > t_rel * norm0:
        basis.append(candidate / norm)
        return True
    return False


def _target_rank(coords: np.ndarray, graph) -> int:
    components = _connected_components(graph)
    rank = 0
    for component in components:
        rank += 3 * len(component) - (5 if _is_linear(coords[list(component)]) else 6)
    return rank


def _connected_components(graph) -> list[tuple[int, ...]]:
    seen: set[int] = set()
    components: list[tuple[int, ...]] = []
    for start in range(graph.natoms):
        if start in seen:
            continue
        stack = [start]
        component = []
        seen.add(start)
        while stack:
            atom = stack.pop()
            component.append(atom)
            for neighbor in graph.adjacency[atom]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return components


def _is_linear(coords: np.ndarray) -> bool:
    if coords.shape[0] <= 2:
        return True
    centered = coords - coords.mean(axis=0)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    return bool(singular_values[1] <= max(1.0e-8, 1.0e-6 * singular_values[0]))


def _primitive_basis(coordinates: Iterable[GICForgePythonCoordinate]) -> tuple[Primitive, ...]:
    basis: list[Primitive] = []
    seen: set[Primitive] = set()
    for coordinate in coordinates:
        for _coefficient, primitive in coordinate.terms:
            if primitive in seen:
                continue
            seen.add(primitive)
            basis.append(primitive)
    return tuple(basis)


def _format_readgic(name: str, terms: tuple[tuple[float, Primitive], ...]) -> str:
    return f"{name}={_format_terms(terms)}"


def _format_terms(terms: tuple[tuple[float, Primitive], ...]) -> str:
    chunks = []
    for coefficient, primitive in terms:
        atom_text = ",".join(str(atom + 1) for atom in primitive.atoms)
        symbol = {
            "bond": "R",
            "angle": "A",
            "linear_bend": "L",
            "dihedral": "D",
            "out_of_plane": "U",
        }[primitive.kind]
        if primitive.kind == "linear_bend":
            atom_text = f"{atom_text},{primitive.mode}"
        chunks.append(f"{coefficient:.10g}*{symbol}({atom_text})")
    return "+".join(chunks)
