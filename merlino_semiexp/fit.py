from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
from io import StringIO

import numpy as np

from geometry.rotational import rotational_constants_MHz
from geometry.structure import Structure
from merlino_core import ScientificValidationError, build_run_manifest
from merlino_fit.survibfit.modify_geom import read_xyz, write_xyz
from merlino_fit.survibfit.pipeline import b_matrix_analytic, build_topology, primitives_from_topology
from merlino_fit.survibfit.primitives import eval_primitives
from merlino_fit.survibfit.transforms import build_u
from merlino_vpt2_vci.internal_gf import gic_labels_from_u, primitive_label

from .contracts import IsotopologueObservation, SemiexperimentalFitRequest


@dataclass(frozen=True)
class SemiexperimentalParameter:
    name: str
    value: float
    sigma: float
    active: bool


@dataclass(frozen=True)
class SemiexperimentalResidual:
    isotopologue: str
    constant: str
    observed_equilibrium_MHz: float
    calculated_MHz: float
    residual_MHz: float


@dataclass(frozen=True)
class SemiexperimentalFitDiagnostics:
    convergence_reason: str
    objective: float
    weighted_rms: float
    reduced_chi_square: float
    rank: int
    condition_number: float
    damping: float
    accepted_steps: int
    rejected_steps: int


@dataclass(frozen=True)
class SemiexperimentalFitResult:
    atoms: tuple[str, ...]
    initial_coordinates_angstrom: np.ndarray
    final_coordinates_angstrom: np.ndarray
    parameters: tuple[SemiexperimentalParameter, ...]
    residuals: tuple[SemiexperimentalResidual, ...]
    covariance: np.ndarray
    correlation: np.ndarray
    jacobian: np.ndarray
    hessian: np.ndarray
    hessian_eigenvalues: np.ndarray
    stationary_point: str
    gic_labels: tuple[str, ...]
    b_matrix: np.ndarray
    iterations: int
    rms_MHz: float
    diagnostics: SemiexperimentalFitDiagnostics
    manifest: Path | None = None


def fit_semiexperimental_geometry(
    request: SemiexperimentalFitRequest,
    *,
    max_iter: int = 12,
    step: float = 1.0e-4,
    damping: float = 1.0e-8,
    max_step: float = 0.25,
    tolerance_MHz: float = 1.0e-6,
    gradient_tolerance: float = 1.0e-8,
    outdir: Path | None = None,
) -> SemiexperimentalFitResult:
    """Fit equilibrium geometry to semiexperimental rotational constants.

    Geometry parameters are Merlino non-redundant GICs. Fixed parameters are
    matched by substring against generated GIC names, or by exact `GICnnn`.
    """
    request.validate()
    atoms, coords0, _comment = read_xyz(Path(request.initial_geometry))
    coords = np.asarray(coords0, dtype=float)
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    _validate_observations(request.observations, len(atoms))

    prims, u_matrix, labels = _gic_model(coords, z_numbers)
    active_mask = _active_mask(labels, request.fixed_parameters)
    if not np.any(active_mask):
        raise ScientificValidationError("All semiexperimental GIC parameters are fixed")

    current_damping = float(damping)
    accepted_steps = 0
    rejected_steps = 0
    convergence_reason = "max_iter"
    previous_objective = None
    for iteration in range(1, max_iter + 1):
        prims, u_matrix, labels = _gic_model(coords, z_numbers)
        active_mask = _active_mask(labels, request.fixed_parameters)
        q = _gic_values(prims, u_matrix, coords)
        calc = _constants_vector(atoms, coords, request.observations)
        obs = _observed_vector(request.observations)
        weights = _weights_vector(request.observations)
        sqrt_weights = np.sqrt(weights)
        residual = obs - calc
        weighted_residual = residual * sqrt_weights
        objective = _objective(weighted_residual)
        jac = _jacobian_constants_wrt_gics(
            atoms,
            coords,
            request.observations,
            prims,
            u_matrix,
            active_mask,
            step=step,
        )
        if np.sqrt(np.mean(residual * residual)) < tolerance_MHz:
            convergence_reason = "rms_tolerance"
            break
        jac_weighted = jac * sqrt_weights[:, None]
        gradient = jac_weighted.T @ weighted_residual
        if float(np.linalg.norm(gradient, ord=np.inf)) < gradient_tolerance:
            convergence_reason = "gradient_tolerance"
            break
        lhs = jac_weighted.T @ jac_weighted + current_damping * np.eye(jac.shape[1])
        dq_active = np.linalg.solve(lhs, gradient)
        dq_active = _limit_step(dq_active, max_step)
        dq = np.zeros_like(q)
        dq[np.where(active_mask)[0]] = dq_active
        candidate, candidate_objective = _line_search_update(
            atoms,
            coords,
            z_numbers,
            request.observations,
            prims,
            u_matrix,
            dq,
            current_objective=objective,
        )
        if candidate_objective < objective:
            coords = candidate
            accepted_steps += 1
            current_damping = max(current_damping / 3.0, 1.0e-14)
            if previous_objective is not None and abs(previous_objective - candidate_objective) < tolerance_MHz * tolerance_MHz:
                convergence_reason = "objective_tolerance"
                break
            previous_objective = candidate_objective
        else:
            rejected_steps += 1
            current_damping = min(current_damping * 10.0, 1.0e12)
    else:
        iteration = max_iter

    prims, u_matrix, labels = _gic_model(coords, z_numbers)
    active_mask = _active_mask(labels, request.fixed_parameters)
    q_final = _gic_values(prims, u_matrix, coords)
    bq = u_matrix.T @ b_matrix_analytic(prims, coords)
    calc = _constants_vector(atoms, coords, request.observations)
    obs = _observed_vector(request.observations)
    residual = obs - calc
    jac = _jacobian_constants_wrt_gics(
        atoms, coords, request.observations, prims, u_matrix, active_mask, step=step
    )
    sqrt_weights = np.sqrt(_weights_vector(request.observations))
    weighted_jac = jac * sqrt_weights[:, None]
    weighted_residual = residual * sqrt_weights
    hessian = _least_squares_hessian(weighted_jac)
    covariance = _covariance(weighted_jac, weighted_residual)
    correlation = _correlation(covariance)
    hessian_eigenvalues = np.linalg.eigvalsh(hessian) if hessian.size else np.array(())
    stationary_point = _stationary_point_type(hessian_eigenvalues)
    diagnostics = _diagnostics(
        weighted_jac,
        weighted_residual,
        convergence_reason=convergence_reason,
        damping=current_damping,
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
    )
    sigmas_active = np.sqrt(np.clip(np.diag(covariance), 0.0, None)) if covariance.size else np.array(())
    parameters = _parameters(labels, q_final, active_mask, sigmas_active)
    residual_rows = _residual_rows(request.observations, calc, obs)
    rms = float(np.sqrt(np.mean(residual * residual))) if residual.size else 0.0
    manifest = None
    if outdir is not None:
        manifest = write_semiexperimental_outputs(
            Path(outdir),
            request,
            atoms,
            coords,
            parameters,
            residual_rows,
            covariance=covariance,
            correlation=correlation,
            hessian=hessian,
            hessian_eigenvalues=hessian_eigenvalues,
            stationary_point=stationary_point,
            diagnostics=diagnostics,
        )
    return SemiexperimentalFitResult(
        atoms=tuple(atoms),
        initial_coordinates_angstrom=np.asarray(coords0, dtype=float),
        final_coordinates_angstrom=coords,
        parameters=parameters,
        residuals=residual_rows,
        covariance=covariance,
        correlation=correlation,
        jacobian=jac,
        hessian=hessian,
        hessian_eigenvalues=hessian_eigenvalues,
        stationary_point=stationary_point,
        gic_labels=labels,
        b_matrix=bq,
        iterations=iteration,
        rms_MHz=rms,
        diagnostics=diagnostics,
        manifest=manifest,
    )


def write_semiexperimental_outputs(
    outdir: Path,
    request: SemiexperimentalFitRequest,
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    parameters: tuple[SemiexperimentalParameter, ...],
    residuals: tuple[SemiexperimentalResidual, ...],
    covariance: np.ndarray | None = None,
    correlation: np.ndarray | None = None,
    hessian: np.ndarray | None = None,
    hessian_eigenvalues: np.ndarray | None = None,
    stationary_point: str = "not_checked",
    diagnostics: SemiexperimentalFitDiagnostics | None = None,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    xyz = outdir / "semiexp_geometry.xyz"
    params = outdir / "semiexp_parameters.csv"
    residual_csv = outdir / "semiexp_residuals.csv"
    covariance_csv = outdir / "semiexp_covariance.csv"
    correlation_csv = outdir / "semiexp_correlation.csv"
    hessian_csv = outdir / "semiexp_hessian.csv"
    hessian_eigs_csv = outdir / "semiexp_hessian_eigenvalues.csv"
    diagnostics_csv = outdir / "semiexp_diagnostics.csv"
    active_names = tuple(p.name for p in parameters if p.active)
    write_xyz(xyz, atoms, coords, comment="Merlino semiexperimental equilibrium geometry")
    params.write_text(parameters_csv(parameters), encoding="utf-8")
    residual_csv.write_text(residuals_csv(residuals), encoding="utf-8")
    covariance_csv.write_text(_matrix_csv(active_names, covariance), encoding="utf-8")
    correlation_csv.write_text(_matrix_csv(active_names, correlation), encoding="utf-8")
    hessian_csv.write_text(_matrix_csv(active_names, hessian), encoding="utf-8")
    hessian_eigs_csv.write_text(_eigenvalues_csv(hessian_eigenvalues), encoding="utf-8")
    diagnostics_csv.write_text(_diagnostics_csv(diagnostics), encoding="utf-8")
    manifest = build_run_manifest(
        workflow="semiexperimental_geometry",
        status="completed",
        run_dir=outdir,
        inputs={"initial_geometry": request.initial_geometry},
        outputs={
            "geometry": xyz,
            "parameters": params,
            "residuals": residual_csv,
            "covariance": covariance_csv,
            "correlation": correlation_csv,
            "hessian": hessian_csv,
            "hessian_eigenvalues": hessian_eigs_csv,
            "diagnostics": diagnostics_csv,
        },
        parameters={
            "fixed_parameters": request.fixed_parameters,
            "stationary_point": stationary_point,
            "convergence_reason": diagnostics.convergence_reason if diagnostics else "not_reported",
        },
        backend={"solver": "python", "coordinate_model": "merlino-gic", "b_matrix": "analytic"},
    )
    return manifest.write(outdir / "semiexp_manifest.json")


def parameters_csv(parameters: tuple[SemiexperimentalParameter, ...]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["name", "value", "sigma", "active"])
    for p in parameters:
        writer.writerow([p.name, f"{p.value:.12g}", f"{p.sigma:.12g}", int(p.active)])
    return stream.getvalue()


def residuals_csv(residuals: tuple[SemiexperimentalResidual, ...]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["isotopologue", "constant", "observed_equilibrium_MHz", "calculated_MHz", "residual_MHz"])
    for r in residuals:
        writer.writerow([
            r.isotopologue,
            r.constant,
            f"{r.observed_equilibrium_MHz:.12g}",
            f"{r.calculated_MHz:.12g}",
            f"{r.residual_MHz:.12g}",
        ])
    return stream.getvalue()


def _matrix_csv(labels: tuple[str, ...], matrix: np.ndarray | None) -> str:
    mat = np.asarray(matrix if matrix is not None else np.zeros((0, 0)), dtype=float)
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["parameter", *labels])
    for label, row in zip(labels, mat):
        writer.writerow([label, *[f"{value:.12g}" for value in row]])
    return stream.getvalue()


def _eigenvalues_csv(values: np.ndarray | None) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["index", "eigenvalue"])
    for idx, value in enumerate(np.asarray(values if values is not None else (), dtype=float), start=1):
        writer.writerow([idx, f"{value:.12g}"])
    return stream.getvalue()


def _diagnostics_csv(diagnostics: SemiexperimentalFitDiagnostics | None) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["key", "value"])
    if diagnostics is None:
        writer.writerow(["status", "not_reported"])
        return stream.getvalue()
    for key, value in diagnostics.__dict__.items():
        writer.writerow([key, value])
    return stream.getvalue()


def _gic_model(coords: np.ndarray, z_numbers: np.ndarray):
    prims = primitives_from_topology(coords, z_numbers, np.deg2rad(170.0), coords_units="angstrom")
    _, _, ringset = build_topology(coords, z_numbers, coords_units="angstrom")
    u_matrix = build_u(prims, coords, Z=z_numbers, ringset=ringset, tol=1.0e-8, fd_step=1.0e-4)
    if u_matrix.size == 0:
        raise ScientificValidationError("Merlino did not generate non-redundant GICs")
    primitive_labels = tuple(primitive_label(p) for p in prims)
    labels = gic_labels_from_u(u_matrix, primitive_labels, threshold=0.20)
    labels = tuple(f"GIC{idx + 1:03d} {label}" for idx, label in enumerate(labels))
    return prims, u_matrix, labels


def _gic_values(prims: object, u_matrix: np.ndarray, coords: np.ndarray) -> np.ndarray:
    return u_matrix.T @ eval_primitives(prims, coords)


def _active_mask(labels: tuple[str, ...], fixed: tuple[str, ...]) -> np.ndarray:
    if not fixed:
        return np.ones(len(labels), dtype=bool)
    mask = []
    fixed_l = tuple(item.lower() for item in fixed)
    for label in labels:
        low = label.lower()
        mask.append(not any(item and item in low for item in fixed_l))
    return np.array(mask, dtype=bool)


def _jacobian_constants_wrt_gics(
    atoms: list[str],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
    prims: object,
    u_matrix: np.ndarray,
    active_mask: np.ndarray,
    *,
    step: float,
) -> np.ndarray:
    active_indices = np.where(active_mask)[0]
    base_q = _gic_values(prims, u_matrix, coords)
    jac = np.zeros((3 * len(observations), len(active_indices)), dtype=float)
    for col, idx in enumerate(active_indices):
        dq = np.zeros_like(base_q)
        dq[idx] = step
        plus = _displace_along_gics(coords, prims, u_matrix, dq)
        dq[idx] = -step
        minus = _displace_along_gics(coords, prims, u_matrix, dq)
        jac[:, col] = (_constants_vector(atoms, plus, observations) - _constants_vector(atoms, minus, observations)) / (2.0 * step)
    return jac


def _displace_along_gics(coords: np.ndarray, prims: object, u_matrix: np.ndarray, dq: np.ndarray) -> np.ndarray:
    bq = u_matrix.T @ b_matrix_analytic(prims, coords)
    dx = np.linalg.pinv(bq, rcond=1.0e-8) @ dq
    return coords + dx.reshape(coords.shape)


def _line_search_update(
    atoms: list[str],
    coords: np.ndarray,
    z_numbers: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
    prims: object,
    u_matrix: np.ndarray,
    dq: np.ndarray,
    *,
    current_objective: float,
) -> tuple[np.ndarray, float]:
    observed = _observed_vector(observations)
    sqrt_weights = np.sqrt(_weights_vector(observations))
    best_coords = coords
    best_objective = current_objective
    for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625):
        candidate = _displace_along_gics(coords, prims, u_matrix, scale * dq)
        try:
            _gic_model(candidate, z_numbers)
        except Exception:
            continue
        residual = (observed - _constants_vector(atoms, candidate, observations)) * sqrt_weights
        objective = _objective(residual)
        if objective < best_objective:
            best_coords = candidate
            best_objective = objective
            break
    return best_coords, best_objective


def _constants_vector(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
) -> np.ndarray:
    values: list[float] = []
    for obs in observations:
        isotopes = _isotopes_for_observation(atoms, obs)
        structure = Structure.from_atoms_coords(list(atoms), [tuple(row) for row in coords], isotopes=isotopes)
        values.extend(rotational_constants_MHz(structure, isotopic=True))
    return np.array(values, dtype=float)


def _observed_vector(observations: tuple[IsotopologueObservation, ...]) -> np.ndarray:
    values: list[float] = []
    for obs in observations:
        values.extend(obs.corrected.as_tuple())
    return np.array(values, dtype=float)


def _weights_vector(observations: tuple[IsotopologueObservation, ...]) -> np.ndarray:
    values: list[float] = []
    for obs in observations:
        values.extend(obs.weights.as_tuple() if obs.weights is not None else (1.0, 1.0, 1.0))
    return np.array(values, dtype=float)


def _residual_rows(
    observations: tuple[IsotopologueObservation, ...],
    calculated: np.ndarray,
    observed: np.ndarray,
) -> tuple[SemiexperimentalResidual, ...]:
    rows = []
    labels = ("A", "B", "C")
    for iso_idx, obs in enumerate(observations):
        for comp_idx, label in enumerate(labels):
            idx = 3 * iso_idx + comp_idx
            rows.append(SemiexperimentalResidual(obs.label, label, float(observed[idx]), float(calculated[idx]), float(observed[idx] - calculated[idx])))
    return tuple(rows)


def _parameters(
    labels: tuple[str, ...],
    values: np.ndarray,
    active_mask: np.ndarray,
    sigmas_active: np.ndarray,
) -> tuple[SemiexperimentalParameter, ...]:
    params = []
    active_counter = 0
    for idx, label in enumerate(labels):
        active = bool(active_mask[idx])
        sigma = float(sigmas_active[active_counter]) if active and active_counter < len(sigmas_active) else 0.0
        if active:
            active_counter += 1
        params.append(SemiexperimentalParameter(label, float(values[idx]), sigma, active))
    return tuple(params)


def _covariance(jac: np.ndarray, residual: np.ndarray) -> np.ndarray:
    if jac.size == 0:
        return np.zeros((0, 0), dtype=float)
    dof = max(jac.shape[0] - jac.shape[1], 1)
    sigma2 = float(residual @ residual) / dof
    return sigma2 * np.linalg.pinv(jac.T @ jac, rcond=1.0e-10)


def _objective(weighted_residual: np.ndarray) -> float:
    return 0.5 * float(weighted_residual @ weighted_residual)


def _limit_step(step: np.ndarray, max_norm: float) -> np.ndarray:
    norm = float(np.linalg.norm(step))
    if max_norm <= 0.0 or norm <= max_norm:
        return step
    return step * (max_norm / norm)


def _diagnostics(
    weighted_jac: np.ndarray,
    weighted_residual: np.ndarray,
    *,
    convergence_reason: str,
    damping: float,
    accepted_steps: int,
    rejected_steps: int,
) -> SemiexperimentalFitDiagnostics:
    if weighted_jac.size:
        singular = np.linalg.svd(weighted_jac, compute_uv=False)
        threshold = max(weighted_jac.shape) * np.finfo(float).eps * (float(singular[0]) if singular.size else 0.0)
        rank = int(np.sum(singular > threshold))
        if singular.size and singular[-1] > threshold:
            condition = float(singular[0] / singular[-1])
        else:
            condition = float("inf")
    else:
        rank = 0
        condition = float("inf")
    objective = _objective(weighted_residual)
    dof = max(weighted_residual.size - weighted_jac.shape[1], 1) if weighted_jac.ndim == 2 else 1
    return SemiexperimentalFitDiagnostics(
        convergence_reason=convergence_reason,
        objective=objective,
        weighted_rms=float(np.sqrt(np.mean(weighted_residual * weighted_residual))) if weighted_residual.size else 0.0,
        reduced_chi_square=float((weighted_residual @ weighted_residual) / dof) if weighted_residual.size else 0.0,
        rank=rank,
        condition_number=condition,
        damping=float(damping),
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
    )


def _least_squares_hessian(weighted_jac: np.ndarray) -> np.ndarray:
    if weighted_jac.size == 0:
        return np.zeros((0, 0), dtype=float)
    return 2.0 * (weighted_jac.T @ weighted_jac)


def _correlation(covariance: np.ndarray) -> np.ndarray:
    if covariance.size == 0:
        return np.zeros((0, 0), dtype=float)
    diag = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    denom = np.outer(diag, diag)
    corr = np.zeros_like(covariance)
    np.divide(covariance, denom, out=corr, where=denom > 0.0)
    return np.clip(corr, -1.0, 1.0)


def _stationary_point_type(eigenvalues: np.ndarray) -> str:
    if eigenvalues.size == 0:
        return "not_checked"
    tol = max(1.0e-10, 1.0e-8 * float(np.max(np.abs(eigenvalues))))
    if np.all(eigenvalues > tol):
        return "minimum"
    if np.any(eigenvalues < -tol):
        return "transition_state_or_saddle"
    return "flat_or_rank_deficient"


def _isotopes_for_observation(atoms: list[str] | tuple[str, ...], obs: IsotopologueObservation) -> list[int | None]:
    isotopes: list[int | None] = [None] * len(atoms)
    for atom_index, isotope_a in obs.substitutions.items():
        if atom_index < 1 or atom_index > len(atoms):
            raise ScientificValidationError(f"Isotopologue {obs.label} substitution atom {atom_index} is out of range")
        isotopes[atom_index - 1] = int(isotope_a)
    return isotopes


def _validate_observations(observations: tuple[IsotopologueObservation, ...], natoms: int) -> None:
    for obs in observations:
        for atom_index in obs.substitutions:
            if atom_index < 1 or atom_index > natoms:
                raise ScientificValidationError(f"Isotopologue {obs.label} substitution atom {atom_index} is out of range")
        if any(value <= 0.0 for value in obs.corrected.as_tuple()):
            raise ScientificValidationError(f"Isotopologue {obs.label} has non-positive equilibrium rotational constants")
        if obs.weights is not None and any(value <= 0.0 for value in obs.weights.as_tuple()):
            raise ScientificValidationError(f"Isotopologue {obs.label} has non-positive least-squares weights")


def _atomic_number(symbol: str) -> int:
    from geometry.elements import atomic_number

    z = atomic_number(symbol)
    if z is None:
        raise ScientificValidationError(f"Unknown element symbol {symbol}")
    return int(z)
