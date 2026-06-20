from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from merlino_fit.survibfit.pipeline import b_matrix, build_topology, primitives_from_topology
from merlino_fit.survibfit.transforms import build_u

from .gaussian_qff import hessian_input_from_gaussian_fchk
from .harmonic import solve_wilson_gf
from .models import HessianInput


@dataclass(frozen=True)
class PEDTable:
    """Potential-energy distribution in non-redundant GIC coordinates."""

    values: np.ndarray
    labels: tuple[str, ...]


@dataclass(frozen=True)
class InternalGFResult:
    frequencies_cm: np.ndarray
    force_constants: np.ndarray
    g_matrix: np.ndarray
    b_matrix: np.ndarray
    u_matrix: np.ndarray
    modes_internal: np.ndarray
    ped: PEDTable
    primitive_labels: tuple[str, ...]
    gic_labels: tuple[str, ...]


def primitive_label(primitive: object) -> str:
    kind = getattr(primitive, "kind")
    atoms = tuple(int(i) + 1 for i in getattr(primitive, "atoms"))
    mode = int(getattr(primitive, "mode", 0))
    suffix = f":{mode}" if mode else ""
    return f"{kind}{suffix}({','.join(str(i) for i in atoms)})"


def gic_labels_from_u(u_matrix: np.ndarray, primitive_labels: tuple[str, ...], *, threshold: float = 0.15) -> tuple[str, ...]:
    labels: list[str] = []
    for col in range(u_matrix.shape[1]):
        terms = []
        for row, coeff in enumerate(u_matrix[:, col]):
            if abs(coeff) < threshold:
                continue
            sign = "+" if coeff >= 0.0 else "-"
            terms.append(f"{sign}{abs(coeff):.3f}*{primitive_labels[row]}")
        labels.append(" ".join(terms) if terms else f"GIC{col + 1}")
    return tuple(labels)


def _mass_inverse(masses_amu: np.ndarray) -> np.ndarray:
    weights = np.repeat(1.0 / np.asarray(masses_amu, dtype=float), 3)
    return np.diag(weights)


def _internal_backtransform(bq: np.ndarray, masses_amu: np.ndarray, g_matrix: np.ndarray) -> np.ndarray:
    minv = _mass_inverse(masses_amu)
    return minv @ bq.T @ np.linalg.pinv(g_matrix, rcond=1.0e-10)


def _ped(force_constants: np.ndarray, modes_internal: np.ndarray, eigenvalues: np.ndarray) -> np.ndarray:
    ped = np.zeros((force_constants.shape[0], modes_internal.shape[1]), dtype=float)
    for mode in range(modes_internal.shape[1]):
        lam = eigenvalues[mode]
        if abs(lam) < 1.0e-14:
            continue
        vector = modes_internal[:, mode]
        raw = vector * (force_constants @ vector) / lam
        total = float(np.sum(np.abs(raw)))
        if total > 0.0:
            ped[:, mode] = 100.0 * np.abs(raw) / total
    return ped


def gf_from_cartesian_hessian_and_merlino_gics(
    cartesian_hessian: np.ndarray,
    coordinates_bohr: np.ndarray,
    atomic_numbers: np.ndarray,
    masses_amu: np.ndarray,
    *,
    fd_step: float = 1.0e-4,
    linear_threshold: float = np.deg2rad(170.0),
) -> InternalGFResult:
    """Run Wilson GF from a Cartesian Hessian and Merlino non-redundant GICs."""
    coords = np.asarray(coordinates_bohr, dtype=float)
    z = np.asarray(atomic_numbers, dtype=int)
    masses = np.asarray(masses_amu, dtype=float)
    hessian = np.asarray(cartesian_hessian, dtype=float)
    if coords.shape != (len(masses), 3):
        raise ValueError("Coordinates and masses have inconsistent dimensions")
    if hessian.shape != (3 * len(masses), 3 * len(masses)):
        raise ValueError("Cartesian Hessian has inconsistent dimensions")

    prims = primitives_from_topology(coords, z, linear_threshold)
    _, _, ringset = build_topology(coords, z, coords_units="au")
    b_prim = b_matrix(prims, coords, fd_step)
    u_matrix = build_u(prims, coords, Z=z, ringset=ringset, tol=1.0e-8, fd_step=fd_step)
    if u_matrix.size == 0:
        raise ValueError("Merlino did not generate non-redundant GICs")

    bq = u_matrix.T @ b_prim
    minv = _mass_inverse(masses)
    g_matrix = bq @ minv @ bq.T
    backtransform = _internal_backtransform(bq, masses, g_matrix)
    force_constants = backtransform.T @ hessian @ backtransform
    force_constants = 0.5 * (force_constants + force_constants.T)

    gf = solve_wilson_gf(force_constants, g_matrix, scale_to_cm=True)
    g_eval, g_vec = np.linalg.eigh(0.5 * (g_matrix + g_matrix.T))
    g_inv_half = (g_vec * (1.0 / np.sqrt(np.clip(g_eval, 1.0e-14, None)))) @ g_vec.T
    modes_internal = g_inv_half @ gf.normal_modes
    ped = _ped(force_constants, modes_internal, gf.eigenvalues)
    primitive_labels = tuple(primitive_label(p) for p in prims)
    gic_labels = gic_labels_from_u(u_matrix, primitive_labels)
    return InternalGFResult(
        frequencies_cm=gf.frequencies_cm,
        force_constants=force_constants,
        g_matrix=g_matrix,
        b_matrix=bq,
        u_matrix=u_matrix,
        modes_internal=modes_internal,
        ped=PEDTable(ped, gic_labels),
        primitive_labels=primitive_labels,
        gic_labels=gic_labels,
    )


def gf_from_hessian_input_with_merlino_gics(input_data: HessianInput) -> InternalGFResult:
    """Run GF/PED from the canonical Merlino harmonic input."""
    input_data.validate()
    return gf_from_cartesian_hessian_and_merlino_gics(
        input_data.cartesian_hessian,
        input_data.cartesian_coordinates_bohr,
        input_data.atomic_numbers,
        input_data.masses_amu,
    )


def gf_from_gaussian_fchk_with_merlino_gics(path: Path) -> InternalGFResult:
    """Gaussian adapter: read FCHK, then run the canonical Merlino GF path."""
    return gf_from_hessian_input_with_merlino_gics(hessian_input_from_gaussian_fchk(path))
