from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
from io import StringIO
import os
import re
import tempfile

import numpy as np

from geometry.inertia import principal_moments
from geometry.physical_constants import Phy, get_physical_constants
from geometry.rotational import rotational_constants_MHz
from geometry.structure import Structure
from merlino_core import ScientificValidationError, build_run_manifest
from merlino_gic import GICDefinition, define_gics_from_cartesian, run_gicforge
from merlino_gic.gic_symmetry import SYMM_INERTIA_TOL as GIC_SYMM_INERTIA_TOL
from merlino_gic.gic_symmetry import SYMM_TOL as GIC_SYMM_TOL
from merlino_core.numerics import damped_normal_step, limit_step, objective, rank_condition
from topology.elements import atomic_symbol
from merlino_fit.topology.pipeline import build_topology_objects
from merlino_fit.survibfit.modify_geom import write_xyz
from merlino_fit.survibfit.pipeline import b_matrix_analytic
from merlino_fit.survibfit.primitives import Primitive, build_primitives, eval_primitives
from merlino_fit.survibfit.symmetry_detector import orient_coords, symmetry_elements_from_geometry
from merlino_fit.survibfit.symmetry_global import primitive_permutation

from .contracts import (
    HYDROGEN_PARAMETER_CONSTRAINT,
    IsotopologueObservation,
    ParameterClassConstraint,
    QMParameterPredicate,
    SemiexperimentalFitRequest,
)
from .geometry_input import read_geometry_input
from .kraitchman import KraitchmanComparison, KraitchmanSeedResult, kraitchman_comparison, kraitchman_seed_geometry
from .cartesian_coordinates import CartesianCoordinateModel, cartesian_symmetry_coordinate_model


ROTATIONAL_COMPONENTS = ("A", "B", "C")
MOMENT_COMPONENTS = ("Ia", "Ib", "Ic")
ROTCONST_TO_MOMENT = (
    get_physical_constants()[Phy.PLANCK]
    / (8.0 * np.pi**2 * get_physical_constants()[Phy.TO_KG] * (1.0e-10) ** 2)
    * 1.0e-6
)


@dataclass(frozen=True)
class SemiexperimentalParameter:
    name: str
    value: float
    sigma: float
    active: bool
    parameter_class: str = ""


@dataclass(frozen=True)
class SemiexperimentalResidual:
    isotopologue: str
    constant: str
    observed_equilibrium_MHz: float
    calculated_MHz: float
    residual_MHz: float


@dataclass(frozen=True)
class SemiexperimentalRotationalConstantComparison:
    isotopologue: str
    component: str
    corrected_experimental_MHz: float
    calculated_MHz: float
    difference_MHz: float


@dataclass(frozen=True)
class SemiexperimentalGeometryParameter:
    kind: str
    label: str
    atom_indices: tuple[int, ...]
    atom_symbols: tuple[str, ...]
    value_angstrom: float | None = None
    value_degree: float | None = None
    sigma_angstrom: float | None = None
    sigma_degree: float | None = None


@dataclass(frozen=True)
class SemiexperimentalFitDiagnostics:
    convergence_reason: str
    objective: float
    weighted_rms: float
    reduced_chi_square: float
    rank: int
    incremental_rank: int
    condition_number: float
    damping: float
    accepted_steps: int
    rejected_steps: int
    max_iterations: int
    n_optimized_parameters: int
    observable: str
    components: tuple[str, ...]
    planar: bool
    auto_pruned_parameters: tuple[str, ...] = ()
    prune_condition_target: float = 0.0
    gicforge_calls: int = 0
    coordinate_model_reuse_steps: int = 0
    trust_radius: float = 0.0
    last_trust_ratio: float = 0.0
    last_line_search_scale: float = 0.0
    b_projector_analytic_refreshes: int = 0
    b_projector_secant_updates: int = 0
    b_projector_secant_rejections: int = 0
    last_b_projector_secant_error: float = 0.0
    parameter_scale_min: float = 1.0
    parameter_scale_max: float = 1.0
    coordinate_model: str = "gic"
    solver: str = "adaptive_lm_trust_region"


@dataclass(frozen=True)
class MeasurementModel:
    observable: str
    components: tuple[str, ...]
    labels: tuple[tuple[str, str], ...]
    observed: np.ndarray
    weights: np.ndarray
    planar: bool


@dataclass(frozen=True)
class LineSearchResult:
    coords: np.ndarray
    q_values: np.ndarray
    objective: float
    accepted: bool
    actual_reduction: float
    predicted_reduction: float
    ratio: float
    scale: float


@dataclass(frozen=True)
class GICProjectorState:
    coords: np.ndarray
    q_values: np.ndarray
    cartesian_from_q: np.ndarray


@dataclass(frozen=True)
class SecantProjectorUpdate:
    cartesian_from_q: np.ndarray | None
    relative_error: float
    accepted: bool


@dataclass
class GICForgeSEBackend:
    atoms: tuple[str, ...]
    root: Path
    counter: int = 0
    last_workdir: Path | None = None
    point_group: str | None = None
    definition: GICDefinition | None = None

    def model(self, coords: np.ndarray):
        if self.definition is not None:
            return self.definition.model()
        self.counter += 1
        workdir = self.root / f"iter_{self.counter:04d}"
        workdir.mkdir(parents=True, exist_ok=True)
        definition = define_gics_from_cartesian(
            self.atoms,
            coords,
            workdir=workdir,
            runner=run_gicforge,
        )
        point_group = _gicforge_point_group(workdir / "provout")
        if self.point_group is None:
            self.point_group = point_group
        elif point_group != self.point_group:
            if _is_symmetry_refinement(self.point_group, point_group):
                self.point_group = point_group
            else:
                raise ScientificValidationError(
                    f"GICForge point group changed from {self.point_group} to {point_group} in {workdir}"
                )
        self.last_workdir = workdir
        self.definition = definition
        return definition.model()


@dataclass(frozen=True)
class SemiexperimentalFitResult:
    atoms: tuple[str, ...]
    initial_coordinates_angstrom: np.ndarray
    final_coordinates_angstrom: np.ndarray
    parameters: tuple[SemiexperimentalParameter, ...]
    geometry_parameters: tuple[SemiexperimentalGeometryParameter, ...]
    residuals: tuple[SemiexperimentalResidual, ...]
    rotational_constants: tuple[SemiexperimentalRotationalConstantComparison, ...]
    kraitchman: tuple[KraitchmanComparison, ...]
    kraitchman_seed: KraitchmanSeedResult | None
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
    max_iter: int | None = None,
    step: float = 1.0e-4,
    damping: float = 1.0e-8,
    max_step: float = 0.25,
    prune_condition: float = 0.0,
    tolerance_MHz: float = 1.0e-6,
    gradient_tolerance: float = 1.0e-8,
    outdir: Path | None = None,
) -> SemiexperimentalFitResult:
    """Fit equilibrium geometry to semiexperimental rotational constants.

    The default working coordinates are Merlino non-redundant GICs. As an
    alternative, `coordinate_model="cartesian_symmetry"` uses a Hessian-free
    translation/rotation-free symmetry-adapted Cartesian displacement basis.
    """
    request.validate()
    if request.coordinate_model == "cartesian_symmetry":
        return _fit_semiexperimental_geometry_cartesian_symmetry(
            request,
            max_iter=max_iter,
            step=step,
            damping=damping,
            max_step=max_step,
            prune_condition=prune_condition,
            tolerance_MHz=tolerance_MHz,
            gradient_tolerance=gradient_tolerance,
            outdir=outdir,
        )
    geometry_input = read_geometry_input(Path(request.initial_geometry))
    atoms = list(geometry_input.atoms)
    coords = np.asarray(geometry_input.coordinates_angstrom, dtype=float)
    coords0 = coords.copy()
    fixed_parameters = _combined_fixed_parameters(request.fixed_parameters, geometry_input.fixed_parameters)
    fixed_gic_patterns = _gic_fixed_patterns(fixed_parameters)
    fixed_primitives = _fixed_primitives_from_patterns(fixed_parameters)
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    _validate_observations(request.observations, len(atoms))
    gicforge_backend = _make_gicforge_backend(tuple(atoms), outdir)

    prims, u_matrix, labels = _gic_model(coords, z_numbers, request, gicforge_backend)
    fixed_primitives = _merge_primitives(
        fixed_primitives,
        _hydrogen_fixed_primitives(atoms, prims, fixed_parameters, coords=coords),
    )
    fixed_primitives = _symmetry_expanded_fixed_primitives(atoms, coords, prims, fixed_primitives)
    reference_gic_signature = _gic_model_signature(labels)
    measurement_model = _build_measurement_model(request, atoms, coords, prims, u_matrix, labels)
    active_mask = _active_mask(labels, fixed_gic_patterns, request.parameter_classes) & _gicforge_a1_mask(labels)
    initial_transform, _initial_names, _initial_classes = _parameter_class_transform(
        labels, active_mask, request.parameter_classes
    )
    initial_transform, _initial_names = _primitive_constrained_transform(
        coords, prims, u_matrix, active_mask, initial_transform, _initial_names, fixed_primitives
    )
    auto_pruned_patterns: tuple[str, ...] = ()
    if prune_condition > 0.0 and initial_transform.shape[1] > 1:
        try:
            initial_jac_gic = _jacobian_constants_wrt_gics(
                atoms,
                coords,
                request,
                prims,
                u_matrix,
                active_mask,
                labels,
                measurement_model,
                step=step,
            )
            initial_weighted_jac = (initial_jac_gic @ initial_transform) * np.sqrt(measurement_model.weights)[:, None]
            auto_pruned_patterns = _weak_parameter_patterns(_initial_names, initial_weighted_jac, prune_condition)
            if auto_pruned_patterns:
                active_mask &= _auto_pruned_active_mask(labels, auto_pruned_patterns)
                initial_transform, _initial_names, _initial_classes = _parameter_class_transform(
                    labels, active_mask, request.parameter_classes
                )
                initial_transform, _initial_names = _primitive_constrained_transform(
                    coords, prims, u_matrix, active_mask, initial_transform, _initial_names, fixed_primitives
                )
        except Exception:
            # Pruning is an observability refinement; unsupported mock/legacy primitives must not block the fit.
            auto_pruned_patterns = ()
    n_optimized_parameters = initial_transform.shape[1]
    loop_max_iter = _resolve_max_iterations(max_iter, n_optimized_parameters) if n_optimized_parameters else 0

    current_damping = max(float(damping), 0.0)
    trust_radius = float(max_step) if max_step > 0.0 else 0.0
    accepted_steps = 0
    rejected_steps = 0
    stalled_rejections = 0
    model_age = 0
    coordinate_model_reuse_steps = 0
    q_initial = _gic_values(prims, u_matrix, coords)
    projector_state = _gic_projector_state(prims, u_matrix, coords, q_initial)
    b_projector_analytic_refreshes = 1
    b_projector_secant_updates = 0
    b_projector_secant_rejections = 0
    last_b_projector_secant_error = 0.0
    parameter_scale_min = 1.0
    parameter_scale_max = 1.0
    last_trust_ratio = 0.0
    last_line_search_scale = 0.0
    convergence_reason = "max_iter" if loop_max_iter else "no_active_totally_symmetric_parameters"
    previous_objective = None
    iteration = 0
    for iteration in range(1, loop_max_iter + 1):
        active_mask = _active_mask(labels, fixed_gic_patterns, request.parameter_classes)
        active_mask &= _gicforge_a1_mask(labels)
        active_mask &= _auto_pruned_active_mask(labels, auto_pruned_patterns)
        q = _gic_values(prims, u_matrix, coords)
        calc = _measurement_vector(atoms, coords, request, q, labels, measurement_model)
        obs = measurement_model.observed
        weights = measurement_model.weights
        sqrt_weights = np.sqrt(weights)
        residual = obs - calc
        weighted_residual = residual * sqrt_weights
        current_objective = objective(weighted_residual)
        jac_gic = _jacobian_constants_wrt_gics(
            atoms,
            coords,
            request,
            prims,
            u_matrix,
            active_mask,
            labels,
            measurement_model,
            step=step,
            cartesian_from_q=projector_state.cartesian_from_q,
        )
        transform, _reduced_names, _class_by_gic = _parameter_class_transform(labels, active_mask, request.parameter_classes)
        transform, _reduced_names = _primitive_constrained_transform(
            coords,
            prims,
            u_matrix,
            active_mask,
            transform,
            _reduced_names,
            fixed_primitives,
            cartesian_from_q=projector_state.cartesian_from_q,
        )
        jac = jac_gic @ transform
        reduced_scales = _reduced_parameter_scales(labels, active_mask, transform)
        if reduced_scales.size:
            parameter_scale_min = min(parameter_scale_min, float(np.min(reduced_scales)))
            parameter_scale_max = max(parameter_scale_max, float(np.max(reduced_scales)))
        if np.sqrt(np.mean(residual * residual)) < tolerance_MHz:
            convergence_reason = "rms_tolerance"
            break
        jac_weighted = jac * sqrt_weights[:, None]
        jac_weighted_scaled = jac_weighted * reduced_scales[None, :] if reduced_scales.size else jac_weighted
        gradient = jac_weighted_scaled.T @ weighted_residual
        if float(np.linalg.norm(gradient, ord=np.inf)) < gradient_tolerance:
            convergence_reason = "gradient_tolerance"
            break
        dq_scaled = _adaptive_lm_step(jac_weighted_scaled, weighted_residual, current_damping, trust_radius)
        dq_reduced = reduced_scales * dq_scaled if reduced_scales.size else dq_scaled
        dq_active = transform @ dq_reduced
        dq = np.zeros_like(q)
        dq[np.where(active_mask)[0]] = dq_active
        line_search = _line_search_update(
            atoms,
            coords,
            request,
            labels,
            measurement_model,
            prims,
            u_matrix,
            dq,
            current_objective=current_objective,
            base_q=q,
            cartesian_from_q=projector_state.cartesian_from_q,
            weighted_residual=weighted_residual,
            jac_weighted=jac_weighted_scaled,
            reduced_step=dq_scaled,
        )
        last_trust_ratio = line_search.ratio
        last_line_search_scale = line_search.scale
        if line_search.accepted:
            previous_coords = coords
            previous_q = q
            next_model_age = model_age + 1
            secant_update = _secant_projector_update(
                projector_state.cartesian_from_q,
                previous_coords,
                previous_q,
                line_search.coords,
                line_search.q_values,
            )
            last_b_projector_secant_error = secant_update.relative_error
            try:
                validation_model = _gic_model(line_search.coords, z_numbers, request, gicforge_backend)
                _validate_gic_model_signature(validation_model[2], reference_gic_signature)
            except Exception:
                rejected_steps += 1
                stalled_rejections += 1
                current_damping, trust_radius = _rejected_trust_update(current_damping, trust_radius, max_step)
                continue
            refresh_required = _should_refresh_gic_model(
                line_search,
                next_model_age,
                secant_relative_error=secant_update.relative_error,
                tolerance_MHz=tolerance_MHz,
                n_observations=len(weighted_residual),
            )
            refreshed_model = None
            if refresh_required:
                refreshed_model = validation_model
            coords = line_search.coords
            accepted_steps += 1
            stalled_rejections = 0
            current_damping, trust_radius = _accepted_trust_update(
                current_damping,
                trust_radius,
                line_search.ratio,
                line_search.scale,
                float(np.linalg.norm(dq_scaled)),
                max_step,
            )
            if previous_objective is not None and abs(previous_objective - line_search.objective) < tolerance_MHz * tolerance_MHz:
                convergence_reason = "objective_tolerance"
                break
            previous_objective = line_search.objective
            if refresh_required and refreshed_model is not None:
                prims, u_matrix, labels = refreshed_model
                refreshed_q = _gic_values(prims, u_matrix, coords)
                projector_state = _gic_projector_state(prims, u_matrix, coords, refreshed_q)
                b_projector_analytic_refreshes += 1
                model_age = 0
            else:
                if not secant_update.accepted or secant_update.cartesian_from_q is None:
                    projector_state = _gic_projector_state(prims, u_matrix, coords, line_search.q_values)
                    b_projector_analytic_refreshes += 1
                    b_projector_secant_rejections += 1
                else:
                    projector_state = GICProjectorState(
                        coords=coords.copy(),
                        q_values=line_search.q_values.copy(),
                        cartesian_from_q=secant_update.cartesian_from_q,
                    )
                    b_projector_secant_updates += 1
                model_age = next_model_age
                coordinate_model_reuse_steps += 1
        else:
            rejected_steps += 1
            stalled_rejections += 1
            current_damping, trust_radius = _rejected_trust_update(current_damping, trust_radius, max_step)
            if model_age:
                prims, u_matrix, labels = _gic_model(coords, z_numbers, request, gicforge_backend)
                _validate_gic_model_signature(labels, reference_gic_signature)
                refreshed_q = _gic_values(prims, u_matrix, coords)
                projector_state = _gic_projector_state(prims, u_matrix, coords, refreshed_q)
                b_projector_analytic_refreshes += 1
                model_age = 0
                stalled_rejections = 0
            if current_damping >= 1.0e12 and stalled_rejections >= 3:
                convergence_reason = "line_search_stalled"
                break
    else:
        iteration = loop_max_iter

    try:
        final_model = _gic_model(coords, z_numbers, request, gicforge_backend)
        _validate_gic_model_signature(final_model[2], reference_gic_signature)
        prims, u_matrix, labels = final_model
    except Exception:
        # Some nearly converged geometries lower only the detected point group at
        # post-processing tolerance.  The last validated GIC model remains the
        # chemically intended coordinate frame for reporting and covariance.
        pass
    active_mask = _active_mask(labels, fixed_gic_patterns, request.parameter_classes)
    active_mask &= _gicforge_a1_mask(labels)
    active_mask &= _auto_pruned_active_mask(labels, auto_pruned_patterns)
    q_final = _gic_values(prims, u_matrix, coords)
    bq = u_matrix.T @ b_matrix_analytic(prims, coords)
    measurement_model = _build_measurement_model(request, atoms, coords, prims, u_matrix, labels)
    calc = _measurement_vector(atoms, coords, request, q_final, labels, measurement_model)
    obs = measurement_model.observed
    residual = obs - calc
    jac_gic = _jacobian_constants_wrt_gics(
        atoms, coords, request, prims, u_matrix, active_mask, labels, measurement_model, step=step
    )
    transform, reduced_names, class_by_gic = _parameter_class_transform(labels, active_mask, request.parameter_classes)
    transform, reduced_names = _primitive_constrained_transform(
        coords, prims, u_matrix, active_mask, transform, reduced_names, fixed_primitives
    )
    jac = jac_gic @ transform
    sqrt_weights = np.sqrt(measurement_model.weights)
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
        max_iterations=loop_max_iter,
        n_optimized_parameters=jac.shape[1],
        observable=measurement_model.observable,
        components=measurement_model.components,
        planar=measurement_model.planar,
        auto_pruned_parameters=auto_pruned_patterns,
        prune_condition_target=prune_condition,
        gicforge_calls=gicforge_backend.counter,
        coordinate_model_reuse_steps=coordinate_model_reuse_steps,
        trust_radius=trust_radius,
        last_trust_ratio=last_trust_ratio,
        last_line_search_scale=last_line_search_scale,
        b_projector_analytic_refreshes=b_projector_analytic_refreshes,
        b_projector_secant_updates=b_projector_secant_updates,
        b_projector_secant_rejections=b_projector_secant_rejections,
        last_b_projector_secant_error=last_b_projector_secant_error,
        parameter_scale_min=parameter_scale_min,
        parameter_scale_max=parameter_scale_max,
    )
    class_by_gic = _mark_auto_pruned_classes(labels, class_by_gic, auto_pruned_patterns)
    parameters = _parameters(labels, q_final, active_mask, transform=transform, covariance=covariance, class_by_gic=class_by_gic)
    residual_rows = _residual_rows(measurement_model, calc, obs)
    rotational_constant_rows = _rotational_constant_rows(atoms, coords, request.observations)
    geometry_parameters = _geometry_parameters(
        atoms,
        coords,
        fit_prims=prims,
        fit_u_matrix=u_matrix,
        active_mask=active_mask,
        transform=transform,
        covariance=covariance,
    )
    kraitchman_rows = kraitchman_comparison(atoms, coords, request.observations)
    kraitchman_seed = kraitchman_seed_geometry(atoms, coords, request.observations, kraitchman_rows)
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
            kraitchman_rows,
            rotational_constants=rotational_constant_rows,
            geometry_parameters=geometry_parameters,
            kraitchman_seed=kraitchman_seed,
            input_fixed_parameters=geometry_input.fixed_parameters,
            effective_parameter_names=reduced_names,
            covariance=covariance,
            correlation=correlation,
            hessian=hessian,
            hessian_eigenvalues=hessian_eigenvalues,
            stationary_point=stationary_point,
            diagnostics=diagnostics,
            measurement_model=measurement_model,
            weighted_jacobian=weighted_jac,
            weighted_residual=weighted_residual,
        )
    return SemiexperimentalFitResult(
        atoms=tuple(atoms),
        initial_coordinates_angstrom=np.asarray(coords0, dtype=float),
        final_coordinates_angstrom=coords,
        parameters=parameters,
        geometry_parameters=geometry_parameters,
        residuals=residual_rows,
        rotational_constants=rotational_constant_rows,
        kraitchman=kraitchman_rows,
        kraitchman_seed=kraitchman_seed,
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


def _fit_semiexperimental_geometry_cartesian_symmetry(
    request: SemiexperimentalFitRequest,
    *,
    max_iter: int | None,
    step: float,
    damping: float,
    max_step: float,
    prune_condition: float,
    tolerance_MHz: float,
    gradient_tolerance: float,
    outdir: Path | None,
) -> SemiexperimentalFitResult:
    geometry_input = read_geometry_input(Path(request.initial_geometry))
    atoms = list(geometry_input.atoms)
    coords = np.asarray(geometry_input.coordinates_angstrom, dtype=float)
    coords0 = coords.copy()
    fixed_parameters = _combined_fixed_parameters(request.fixed_parameters, geometry_input.fixed_parameters)
    fixed_mode_patterns = _gic_fixed_patterns(fixed_parameters)
    fixed_primitives = _fixed_primitives_from_patterns(fixed_parameters)
    _validate_observations(request.observations, len(atoms))
    mode_model = cartesian_symmetry_coordinate_model(tuple(atoms), coords0)
    labels = mode_model.labels
    constraint_prims = tuple(_constraint_primitive_pool(atoms, (), coords))
    fixed_primitives = _merge_primitives(
        fixed_primitives,
        _hydrogen_fixed_primitives(atoms, constraint_prims, fixed_parameters, coords=coords),
    )
    fixed_primitives = _symmetry_expanded_fixed_primitives(atoms, coords, constraint_prims, fixed_primitives)

    measurement_model = _build_measurement_model_cartesian_basis(
        request,
        atoms,
        coords,
        labels,
        mode_model.cartesian_from_q,
    )
    active_mask = _active_mask(labels, fixed_mode_patterns, request.parameter_classes)
    active_mask &= mode_model.active_totally_symmetric_mask
    transform, reduced_names, class_by_mode = _parameter_class_transform(labels, active_mask, request.parameter_classes)
    transform, reduced_names = _primitive_constrained_cartesian_transform(
        coords,
        mode_model.cartesian_from_q,
        active_mask,
        transform,
        reduced_names,
        fixed_primitives,
    )
    auto_pruned_patterns: tuple[str, ...] = ()
    if prune_condition > 0.0 and transform.shape[1] > 1:
        jac_modes = _jacobian_constants_wrt_cartesian_basis(
            atoms,
            coords,
            request,
            labels,
            measurement_model,
            mode_model.cartesian_from_q,
        )
        weighted = (_active_coordinate_jacobian(jac_modes, active_mask) @ transform) * np.sqrt(measurement_model.weights)[:, None]
        auto_pruned_patterns = _weak_parameter_patterns(reduced_names, weighted, prune_condition)
        if auto_pruned_patterns:
            active_mask &= _auto_pruned_active_mask(labels, auto_pruned_patterns)
            transform, reduced_names, class_by_mode = _parameter_class_transform(
                labels, active_mask, request.parameter_classes
            )
            transform, reduced_names = _primitive_constrained_cartesian_transform(
                coords,
                mode_model.cartesian_from_q,
                active_mask,
                transform,
                reduced_names,
                fixed_primitives,
            )

    n_optimized_parameters = transform.shape[1]
    loop_max_iter = _resolve_max_iterations(max_iter, n_optimized_parameters) if n_optimized_parameters else 0
    current_damping = max(float(damping), 0.0)
    trust_radius = float(max_step) if max_step > 0.0 else 0.0
    accepted_steps = 0
    rejected_steps = 0
    stalled_rejections = 0
    parameter_scale_min = 1.0
    parameter_scale_max = 1.0
    last_trust_ratio = 0.0
    last_line_search_scale = 0.0
    previous_objective = None
    convergence_reason = "max_iter" if loop_max_iter else "no_active_totally_symmetric_cartesian_coordinates"
    iteration = 0

    for iteration in range(1, loop_max_iter + 1):
        active_mask = _active_mask(labels, fixed_mode_patterns, request.parameter_classes)
        active_mask &= mode_model.active_totally_symmetric_mask
        active_mask &= _auto_pruned_active_mask(labels, auto_pruned_patterns)
        q = mode_model.values(coords)
        calc = _measurement_vector(atoms, coords, request, q, labels, measurement_model)
        obs = measurement_model.observed
        weights = measurement_model.weights
        sqrt_weights = np.sqrt(weights)
        residual = obs - calc
        weighted_residual = residual * sqrt_weights
        current_objective = objective(weighted_residual)
        jac_modes = _jacobian_constants_wrt_cartesian_basis(
            atoms,
            coords,
            request,
            labels,
            measurement_model,
            mode_model.cartesian_from_q,
        )
        transform, reduced_names, class_by_mode = _parameter_class_transform(labels, active_mask, request.parameter_classes)
        transform, reduced_names = _primitive_constrained_cartesian_transform(
            coords,
            mode_model.cartesian_from_q,
            active_mask,
            transform,
            reduced_names,
            fixed_primitives,
        )
        jac = _active_coordinate_jacobian(jac_modes, active_mask) @ transform
        reduced_scales = np.ones(jac.shape[1], dtype=float)
        if np.sqrt(np.mean(residual * residual)) < tolerance_MHz:
            convergence_reason = "rms_tolerance"
            break
        jac_weighted = jac * sqrt_weights[:, None]
        jac_weighted_scaled = jac_weighted * reduced_scales[None, :] if reduced_scales.size else jac_weighted
        gradient = jac_weighted_scaled.T @ weighted_residual
        if float(np.linalg.norm(gradient, ord=np.inf)) < gradient_tolerance:
            convergence_reason = "gradient_tolerance"
            break
        dq_scaled = _adaptive_lm_step(jac_weighted_scaled, weighted_residual, current_damping, trust_radius)
        dq_reduced = reduced_scales * dq_scaled if reduced_scales.size else dq_scaled
        dq_active = transform @ dq_reduced
        dq = np.zeros_like(q)
        dq[np.where(active_mask)[0]] = dq_active
        line_search = _line_search_update_cartesian_basis(
            atoms,
            coords,
            request,
            labels,
            measurement_model,
            mode_model,
            dq,
            current_objective=current_objective,
            base_q=q,
            weighted_residual=weighted_residual,
            jac_weighted=jac_weighted_scaled,
            reduced_step=dq_scaled,
        )
        last_trust_ratio = line_search.ratio
        last_line_search_scale = line_search.scale
        if line_search.accepted:
            coords = line_search.coords
            accepted_steps += 1
            stalled_rejections = 0
            current_damping, trust_radius = _accepted_trust_update(
                current_damping,
                trust_radius,
                line_search.ratio,
                line_search.scale,
                float(np.linalg.norm(dq_scaled)),
                max_step,
            )
            if previous_objective is not None and abs(previous_objective - line_search.objective) < tolerance_MHz * tolerance_MHz:
                convergence_reason = "objective_tolerance"
                break
            previous_objective = line_search.objective
        else:
            rejected_steps += 1
            stalled_rejections += 1
            current_damping, trust_radius = _rejected_trust_update(current_damping, trust_radius, max_step)
            if current_damping >= 1.0e12 and stalled_rejections >= 3:
                convergence_reason = "line_search_stalled"
                break
    else:
        iteration = loop_max_iter

    active_mask = _active_mask(labels, fixed_mode_patterns, request.parameter_classes)
    active_mask &= mode_model.active_totally_symmetric_mask
    active_mask &= _auto_pruned_active_mask(labels, auto_pruned_patterns)
    q_final = mode_model.values(coords)
    measurement_model = _build_measurement_model_cartesian_basis(
        request,
        atoms,
        coords,
        labels,
        mode_model.cartesian_from_q,
    )
    calc = _measurement_vector(atoms, coords, request, q_final, labels, measurement_model)
    obs = measurement_model.observed
    residual = obs - calc
    jac_modes = _jacobian_constants_wrt_cartesian_basis(
        atoms,
        coords,
        request,
        labels,
        measurement_model,
        mode_model.cartesian_from_q,
    )
    transform, reduced_names, class_by_mode = _parameter_class_transform(labels, active_mask, request.parameter_classes)
    transform, reduced_names = _primitive_constrained_cartesian_transform(
        coords,
        mode_model.cartesian_from_q,
        active_mask,
        transform,
        reduced_names,
        fixed_primitives,
    )
    jac = _active_coordinate_jacobian(jac_modes, active_mask) @ transform
    sqrt_weights = np.sqrt(measurement_model.weights)
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
        max_iterations=loop_max_iter,
        n_optimized_parameters=jac.shape[1],
        observable=measurement_model.observable,
        components=measurement_model.components,
        planar=measurement_model.planar,
        auto_pruned_parameters=auto_pruned_patterns,
        prune_condition_target=prune_condition,
        trust_radius=trust_radius,
        last_trust_ratio=last_trust_ratio,
        last_line_search_scale=last_line_search_scale,
        parameter_scale_min=parameter_scale_min,
        parameter_scale_max=parameter_scale_max,
        coordinate_model=request.coordinate_model,
    )
    class_by_mode = _mark_auto_pruned_classes(labels, class_by_mode, auto_pruned_patterns)
    parameters = _parameters(labels, q_final, active_mask, transform=transform, covariance=covariance, class_by_gic=class_by_mode)
    residual_rows = _residual_rows(measurement_model, calc, obs)
    rotational_constant_rows = _rotational_constant_rows(atoms, coords, request.observations)
    cartesian_from_parameters = _cartesian_from_reduced_coordinates(mode_model.cartesian_from_q, active_mask, transform)
    geometry_parameters = _geometry_parameters(
        atoms,
        coords,
        cartesian_from_parameters=cartesian_from_parameters,
        covariance=covariance,
    )
    kraitchman_rows = kraitchman_comparison(atoms, coords, request.observations)
    kraitchman_seed = kraitchman_seed_geometry(atoms, coords, request.observations, kraitchman_rows)
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
            kraitchman_rows,
            rotational_constants=rotational_constant_rows,
            geometry_parameters=geometry_parameters,
            kraitchman_seed=kraitchman_seed,
            input_fixed_parameters=geometry_input.fixed_parameters,
            effective_parameter_names=reduced_names,
            covariance=covariance,
            correlation=correlation,
            hessian=hessian,
            hessian_eigenvalues=hessian_eigenvalues,
            stationary_point=stationary_point,
            diagnostics=diagnostics,
            measurement_model=measurement_model,
            weighted_jacobian=weighted_jac,
            weighted_residual=weighted_residual,
        )
    return SemiexperimentalFitResult(
        atoms=tuple(atoms),
        initial_coordinates_angstrom=np.asarray(coords0, dtype=float),
        final_coordinates_angstrom=coords,
        parameters=parameters,
        geometry_parameters=geometry_parameters,
        residuals=residual_rows,
        rotational_constants=rotational_constant_rows,
        kraitchman=kraitchman_rows,
        kraitchman_seed=kraitchman_seed,
        covariance=covariance,
        correlation=correlation,
        jacobian=jac,
        hessian=hessian,
        hessian_eigenvalues=hessian_eigenvalues,
        stationary_point=stationary_point,
        gic_labels=labels,
        b_matrix=mode_model.cartesian_from_q.T,
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
    kraitchman: tuple[KraitchmanComparison, ...] = (),
    rotational_constants: tuple[SemiexperimentalRotationalConstantComparison, ...] | None = None,
    geometry_parameters: tuple[SemiexperimentalGeometryParameter, ...] | None = None,
    kraitchman_seed: KraitchmanSeedResult | None = None,
    effective_parameter_names: tuple[str, ...] = (),
    covariance: np.ndarray | None = None,
    correlation: np.ndarray | None = None,
    hessian: np.ndarray | None = None,
    hessian_eigenvalues: np.ndarray | None = None,
    stationary_point: str = "not_checked",
    diagnostics: SemiexperimentalFitDiagnostics | None = None,
    input_fixed_parameters: tuple[str, ...] = (),
    measurement_model: MeasurementModel | None = None,
    weighted_jacobian: np.ndarray | None = None,
    weighted_residual: np.ndarray | None = None,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    xyz = outdir / "semiexp_geometry.xyz"
    params = outdir / "semiexp_parameters.csv"
    geometry_params = outdir / "semiexp_geometry_parameters.csv"
    residual_csv = outdir / "semiexp_residuals.csv"
    rotconst_csv = outdir / "semiexp_rotational_constants.csv"
    text_report = outdir / "semiexp_report.txt"
    kraitchman_csv = outdir / "semiexp_kraitchman.csv"
    kraitchman_xyz = outdir / "semiexp_kraitchman_geometry.xyz"
    covariance_csv = outdir / "semiexp_covariance.csv"
    correlation_csv = outdir / "semiexp_correlation.csv"
    hessian_csv = outdir / "semiexp_hessian.csv"
    hessian_eigs_csv = outdir / "semiexp_hessian_eigenvalues.csv"
    diagnostics_csv = outdir / "semiexp_diagnostics.csv"
    influence_csv = outdir / "semiexp_influence.csv"
    high_correlation_csv = outdir / "semiexp_high_correlations.csv"
    active_names = effective_parameter_names or _effective_parameter_names(parameters)
    geometry_rows = geometry_parameters if geometry_parameters is not None else _geometry_parameters(atoms, coords)
    rotconst_rows = (
        rotational_constants
        if rotational_constants is not None
        else _rotational_constant_rows(atoms, np.asarray(coords, dtype=float), request.observations)
    )
    fixed_parameters = _combined_fixed_parameters(request.fixed_parameters, input_fixed_parameters)
    write_xyz(xyz, atoms, coords, comment="Merlino semiexperimental equilibrium geometry")
    params.write_text(parameters_csv(parameters), encoding="utf-8")
    geometry_params.write_text(geometry_parameters_csv(geometry_rows), encoding="utf-8")
    residual_csv.write_text(residuals_csv(residuals), encoding="utf-8")
    rotconst_csv.write_text(rotational_constants_csv(rotconst_rows), encoding="utf-8")
    text_report.write_text(
        semiexperimental_text_report(
            request,
            parameters,
            geometry_rows,
            residuals,
            rotconst_rows,
            diagnostics=diagnostics,
            stationary_point=stationary_point,
            fixed_parameters=fixed_parameters,
        ),
        encoding="utf-8",
    )
    kraitchman_csv.write_text(kraitchman_csv_rows(kraitchman), encoding="utf-8")
    if kraitchman_seed is not None:
        write_xyz(
            kraitchman_xyz,
            atoms,
            kraitchman_seed.coordinates_angstrom,
            comment=f"Merlino Kraitchman substitution geometry; method={kraitchman_seed.method}; principal-axis frame",
        )
    covariance_csv.write_text(_matrix_csv(active_names, covariance), encoding="utf-8")
    correlation_csv.write_text(_matrix_csv(active_names, correlation), encoding="utf-8")
    hessian_csv.write_text(_matrix_csv(active_names, hessian), encoding="utf-8")
    hessian_eigs_csv.write_text(_eigenvalues_csv(hessian_eigenvalues), encoding="utf-8")
    diagnostics_csv.write_text(_diagnostics_csv(diagnostics), encoding="utf-8")
    influence_csv.write_text(
        _influence_csv(measurement_model, residuals, weighted_jacobian, weighted_residual),
        encoding="utf-8",
    )
    high_correlation_csv.write_text(_high_correlations_csv(active_names, correlation), encoding="utf-8")
    manifest_inputs = {"initial_geometry": request.initial_geometry}
    if request.coordinate_model == "cartesian_symmetry":
        coordinate_generation = {
            "primitive_source": "Hessian-free Cartesian displacement basis from the parent geometry",
            "reduction": "ordinary Cartesian translations and rotations projected out",
            "symmetry": "Cartesian displacement basis projected with detected point-group irreps",
            "active_subspace": "totally symmetric symmetry-adapted Cartesian displacements only",
            "ring_coordinates": "not used as working coordinates; final primitive internals are reported from final Cartesian geometry",
            "line_search": "trust-region LM with fixed symmetry-Cartesian basis; no GIC B projector is required",
            "restart_policy": "restart jobs rebuild the symmetry-Cartesian basis from the parent geometry",
        }
        backend_coordinate_model = "symmetry-cartesian"
        b_matrix_description = "not required for working-coordinate updates"
    else:
        coordinate_generation = {
            "primitive_source": "GICForge definition utility run once per fresh fit; frozen GIC schema reused for all iterations",
            "reduction": "primitive stretches plus non-redundant non-stretch GIC transform",
            "symmetry": "GICForge/symm.f point group with deterministic final GIC irrep assignment",
            "active_subspace": "GICForge-assigned A1 coordinates only",
            "ring_coordinates": "GICForge ring deformation and puckering coordinates",
            "gicforge_iterations": str(outdir / "gicforge_iterations"),
            "line_search": "trust-region LM with frozen GIC schema, analytic B rebuilds when required, and secant-updated B projector between rebuilds",
            "restart_policy": "restart jobs rebuild the GIC schema; ordinary iterations reuse the saved schema",
        }
        backend_coordinate_model = "gicforge-frozen-definition"
        b_matrix_description = "analytic from frozen GIC definition"
    outputs = {
        "geometry": xyz,
        "parameters": params,
        "geometry_parameters": geometry_params,
        "residuals": residual_csv,
        "rotational_constants": rotconst_csv,
        "text_report": text_report,
        "kraitchman": kraitchman_csv,
        "covariance": covariance_csv,
        "correlation": correlation_csv,
        "hessian": hessian_csv,
        "hessian_eigenvalues": hessian_eigs_csv,
        "diagnostics": diagnostics_csv,
        "influence": influence_csv,
        "high_correlations": high_correlation_csv,
    }
    if kraitchman_seed is not None:
        outputs["kraitchman_geometry"] = kraitchman_xyz
    manifest = build_run_manifest(
        workflow="semiexperimental_geometry",
        status="completed",
        run_dir=outdir,
        inputs=manifest_inputs,
        outputs=outputs,
        parameters={
            "fixed_parameters": fixed_parameters,
            "input_fixed_parameters": input_fixed_parameters,
            "parameter_classes": tuple(
                {"name": item.name, "patterns": item.patterns, "mode": item.mode}
                for item in request.parameter_classes
            ),
            "stationary_point": stationary_point,
            "convergence_reason": diagnostics.convergence_reason if diagnostics else "not_reported",
            "coordinate_model": diagnostics.coordinate_model if diagnostics else request.coordinate_model,
            "observable": diagnostics.observable if diagnostics else request.observable,
            "rotational_components": diagnostics.components if diagnostics else request.rotational_components,
            "isotopologues": tuple(obs.label for obs in request.observations),
            "n_isotopologues": len(request.observations),
            "n_qm_predicates": len(request.qm_predicates),
            "n_gic_parameters": len(parameters),
            "n_working_parameters": len(parameters),
            "n_effective_parameters": len(active_names),
            "n_active_gic_parameters": sum(1 for item in parameters if item.active),
            "n_active_working_parameters": sum(1 for item in parameters if item.active),
            "auto_pruned_parameters": diagnostics.auto_pruned_parameters if diagnostics else (),
            "prune_condition_target": diagnostics.prune_condition_target if diagnostics else 0.0,
            "max_iterations": diagnostics.max_iterations if diagnostics else None,
            "n_kraitchman_rows": len(kraitchman),
            "kraitchman_seed_method": kraitchman_seed.method if kraitchman_seed else "not_available",
            "n_kraitchman_seed_atoms": len(kraitchman_seed.fitted_atom_indices) if kraitchman_seed else 0,
            "rank": diagnostics.rank if diagnostics else None,
            "incremental_rank": diagnostics.incremental_rank if diagnostics else None,
            "condition_number": diagnostics.condition_number if diagnostics else None,
            "weighted_rms": diagnostics.weighted_rms if diagnostics else None,
            "reduced_chi_square": diagnostics.reduced_chi_square if diagnostics else None,
            "gicforge_calls": diagnostics.gicforge_calls if diagnostics else None,
            "coordinate_model_reuse_steps": diagnostics.coordinate_model_reuse_steps if diagnostics else None,
            "trust_radius": diagnostics.trust_radius if diagnostics else None,
            "last_trust_ratio": diagnostics.last_trust_ratio if diagnostics else None,
            "last_line_search_scale": diagnostics.last_line_search_scale if diagnostics else None,
            "b_projector_analytic_refreshes": diagnostics.b_projector_analytic_refreshes if diagnostics else None,
            "b_projector_secant_updates": diagnostics.b_projector_secant_updates if diagnostics else None,
            "b_projector_secant_rejections": diagnostics.b_projector_secant_rejections if diagnostics else None,
            "last_b_projector_secant_error": diagnostics.last_b_projector_secant_error if diagnostics else None,
            "parameter_scale_min": diagnostics.parameter_scale_min if diagnostics else None,
            "parameter_scale_max": diagnostics.parameter_scale_max if diagnostics else None,
            "coordinate_generation": coordinate_generation,
        },
        backend={
            "solver": "python-orchestrated adaptive trust-region Levenberg-Marquardt",
            "coordinate_model": backend_coordinate_model,
            "b_matrix": b_matrix_description,
            "fortran77_role": "validated numerical kernels only",
            "fortran77_source": "fortran/semiexp/semiexp_core.f",
        },
        messages=[
            "Semiexperimental workflow is orchestrated in Python.",
            "Fortran77 semiexp code is kept as an independent validated numerical-kernel layer.",
        ],
    )
    return manifest.write(outdir / "semiexp_manifest.json")


def parameters_csv(parameters: tuple[SemiexperimentalParameter, ...]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["name", "value", "sigma", "active", "parameter_class"])
    for p in parameters:
        writer.writerow([p.name, f"{p.value:.12g}", f"{p.sigma:.12g}", int(p.active), p.parameter_class])
    return stream.getvalue()


def geometry_parameters_csv(parameters: tuple[SemiexperimentalGeometryParameter, ...]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow([
        "kind",
        "label",
        "atoms",
        "symbols",
        "value_angstrom",
        "sigma_angstrom",
        "value_degree",
        "sigma_degree",
    ])
    for item in parameters:
        writer.writerow([
            item.kind,
            item.label,
            "-".join(str(idx) for idx in item.atom_indices),
            "-".join(item.atom_symbols),
            "" if item.value_angstrom is None else f"{item.value_angstrom:.12g}",
            "" if item.sigma_angstrom is None else f"{item.sigma_angstrom:.12g}",
            "" if item.value_degree is None else f"{item.value_degree:.12g}",
            "" if item.sigma_degree is None else f"{item.sigma_degree:.12g}",
        ])
    return stream.getvalue()


def _effective_parameter_names(parameters: tuple[SemiexperimentalParameter, ...]) -> tuple[str, ...]:
    names: list[str] = []
    seen: set[str] = set()
    for parameter in parameters:
        if not parameter.active:
            continue
        name = parameter.parameter_class or parameter.name
        if name not in seen:
            names.append(name)
            seen.add(name)
    return tuple(names)


def _combined_fixed_parameters(
    explicit_fixed: tuple[str, ...],
    input_fixed: tuple[str, ...] = (),
) -> tuple[str, ...]:
    result: list[str] = []
    seen: set[str] = set()
    for item in (*explicit_fixed, *input_fixed):
        text = str(item).strip()
        if text and text not in seen:
            result.append(text)
            seen.add(text)
    return tuple(result)


def residuals_csv(residuals: tuple[SemiexperimentalResidual, ...]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["isotopologue", "observable", "observed", "calculated", "residual"])
    for r in residuals:
        writer.writerow([
            r.isotopologue,
            r.constant,
            f"{r.observed_equilibrium_MHz:.12g}",
            f"{r.calculated_MHz:.12g}",
            f"{r.residual_MHz:.12g}",
        ])
    return stream.getvalue()


def rotational_constants_csv(rows: tuple[SemiexperimentalRotationalConstantComparison, ...]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow([
        "isotopologue",
        "component",
        "corrected_experimental_MHz",
        "calculated_MHz",
        "difference_MHz",
    ])
    for item in rows:
        writer.writerow([
            item.isotopologue,
            item.component,
            f"{item.corrected_experimental_MHz:.12g}",
            f"{item.calculated_MHz:.12g}",
            f"{item.difference_MHz:.12g}",
        ])
    return stream.getvalue()


def semiexperimental_text_report(
    request: SemiexperimentalFitRequest,
    parameters: tuple[SemiexperimentalParameter, ...],
    geometry_parameters: tuple[SemiexperimentalGeometryParameter, ...],
    residuals: tuple[SemiexperimentalResidual, ...],
    rotational_constants: tuple[SemiexperimentalRotationalConstantComparison, ...],
    *,
    diagnostics: SemiexperimentalFitDiagnostics | None = None,
    stationary_point: str = "not_checked",
    fixed_parameters: tuple[str, ...] = (),
) -> str:
    lines: list[str] = [
        "SEFIT TEXT OUTPUT v1",
        "Merlino semiexperimental equilibrium-geometry fit",
        "=" * 72,
        "",
        "[method]",
        f"program = Merlino SEfit",
        f"method = semiexperimental equilibrium-geometry least squares",
        f"solver = {diagnostics.solver if diagnostics is not None else 'adaptive_lm_trust_region'}",
        f"coordinate_model = {request.coordinate_model}",
        f"coordinate_basis = {_coordinate_model_description(request.coordinate_model)}",
        f"initial_geometry = {request.initial_geometry}",
        f"observable = {diagnostics.observable if diagnostics is not None else request.observable}",
        f"components = {','.join(diagnostics.components) if diagnostics is not None else request.rotational_components}",
        f"isotopologues = {', '.join(obs.label for obs in request.observations)}",
        f"stationary_point = {stationary_point}",
        "",
        "[constraints]",
    ]
    if fixed_parameters:
        lines.extend(f"fixed_parameter = {item}" for item in fixed_parameters)
    else:
        lines.append("fixed_parameter = none")
    if request.parameter_classes:
        for item in request.parameter_classes:
            lines.append(f"parameter_class = {item.name}; mode={item.mode}; patterns={'|'.join(item.patterns)}")
    else:
        lines.append("parameter_class = none")
    if request.qm_predicates:
        for item in request.qm_predicates:
            lines.append(
                f"qm_predicate = {item.label_pattern}; value={item.value:.12g}; "
                f"sigma={item.sigma:.12g}; source={item.source}"
            )
    else:
        lines.append("qm_predicate = none")
    lines.extend(["", "[fit_statistics]"])
    if diagnostics is not None:
        lines.extend(
            [
                f"convergence = {diagnostics.convergence_reason}",
                f"iterations = {diagnostics.accepted_steps + diagnostics.rejected_steps}",
                f"accepted_steps = {diagnostics.accepted_steps}",
                f"rejected_steps = {diagnostics.rejected_steps}",
                f"max_iterations = {diagnostics.max_iterations}",
                f"n_optimized_parameters = {diagnostics.n_optimized_parameters}",
                f"objective = {diagnostics.objective:.12g}",
                f"weighted_rms = {diagnostics.weighted_rms:.12g}",
                f"reduced_chi_square = {diagnostics.reduced_chi_square:.12g}",
                f"rank = {diagnostics.rank}",
                f"incremental_rank = {diagnostics.incremental_rank}",
                f"condition_number = {diagnostics.condition_number:.12g}",
                f"damping = {diagnostics.damping:.12g}",
                f"trust_radius = {diagnostics.trust_radius:.12g}",
                f"last_trust_ratio = {diagnostics.last_trust_ratio:.12g}",
                f"last_line_search_scale = {diagnostics.last_line_search_scale:.12g}",
                f"gicforge_calls = {diagnostics.gicforge_calls}",
                f"coordinate_model_reuse_steps = {diagnostics.coordinate_model_reuse_steps}",
                f"b_projector_analytic_refreshes = {diagnostics.b_projector_analytic_refreshes}",
                f"b_projector_secant_updates = {diagnostics.b_projector_secant_updates}",
                f"b_projector_secant_rejections = {diagnostics.b_projector_secant_rejections}",
                f"last_b_projector_secant_error = {diagnostics.last_b_projector_secant_error:.12g}",
                f"parameter_scale_min = {diagnostics.parameter_scale_min:.12g}",
                f"parameter_scale_max = {diagnostics.parameter_scale_max:.12g}",
            ]
        )
    else:
        lines.append("statistics = not_available")

    lines.extend(["", "[working_coordinates]", f"coordinate_count = {len(parameters)}"])
    lines.append("index active class value sigma label")
    for idx, item in enumerate(parameters, start=1):
        lines.append(
            " ".join(
                (
                    str(idx),
                    "yes" if item.active else "no",
                    item.parameter_class or "-",
                    f"{item.value:.12g}",
                    f"{item.sigma:.12g}",
                    item.name,
                )
            )
        )

    lines.extend(["", "[primitive_internal_coordinates]", "Final topological geometry"])
    lines.append("kind label atoms symbols value sigma unit")
    for item in geometry_parameters:
        atoms = "-".join(str(idx) for idx in item.atom_indices)
        symbols = "-".join(item.atom_symbols)
        if item.value_angstrom is not None:
            value = item.value_angstrom
            sigma = item.sigma_angstrom
            unit = "Angstrom"
        else:
            value = item.value_degree
            sigma = item.sigma_degree
            unit = "degree"
        lines.append(
            " ".join(
                (
                    item.kind,
                    item.label,
                    atoms,
                    symbols,
                    "" if value is None else f"{value:.12g}",
                    "" if sigma is None else f"{sigma:.12g}",
                    unit,
                )
            )
        )

    lines.extend(["", "[rotational_constants]", "Rotational constants (MHz)"])
    lines.append("isotopologue component corrected_experimental_MHz calculated_MHz exp_minus_calc_MHz")
    for item in rotational_constants:
        lines.append(
            " ".join(
                (
                    item.isotopologue,
                    item.component,
                    f"{item.corrected_experimental_MHz:.12g}",
                    f"{item.calculated_MHz:.12g}",
                    f"{item.difference_MHz:.12g}",
                )
            )
        )

    lines.extend(["", "[fit_residuals]"])
    lines.append("isotopologue observable observed calculated residual")
    for item in residuals:
        lines.append(
            " ".join(
                (
                    item.isotopologue,
                    item.constant,
                    f"{item.observed_equilibrium_MHz:.12g}",
                    f"{item.calculated_MHz:.12g}",
                    f"{item.residual_MHz:.12g}",
                )
            )
        )
    return "\n".join(lines) + "\n"


def _coordinate_model_description(coordinate_model: str) -> str:
    if coordinate_model == "cartesian_symmetry":
        return "totally symmetric Hessian-free symmetry-adapted Cartesian displacements"
    return "GICForge non-redundant symmetry-adapted GICs; active subspace is totally symmetric"


def _geometry_parameters(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    *,
    fit_prims: object | None = None,
    fit_u_matrix: np.ndarray | None = None,
    active_mask: np.ndarray | None = None,
    transform: np.ndarray | None = None,
    cartesian_from_parameters: np.ndarray | None = None,
    covariance: np.ndarray | None = None,
) -> tuple[SemiexperimentalGeometryParameter, ...]:
    coords = np.asarray(coords, dtype=float)
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    try:
        _continuous, graph, _ringset, _synthons, _aromaticity = build_topology_objects(coords, z_numbers)
    except Exception as exc:
        raise ScientificValidationError(f"Cannot build final geometry parameter table: {exc}") from exc

    specs: list[tuple[str, str, tuple[int, ...], tuple[str, ...], Primitive, float]] = []
    for i, j in sorted(tuple(sorted(pair)) for pair in graph.bonds):
        label = f"R({i + 1},{j + 1})"
        symbols = (str(atoms[i]), str(atoms[j]))
        specs.append(("bond", label, (i + 1, j + 1), symbols, Primitive("bond", (i, j)), 1.0))

    for center in range(len(atoms)):
        neighbors = sorted(graph.adjacency[center])
        for pos, left in enumerate(neighbors):
            for right in neighbors[pos + 1 :]:
                label = f"A({left + 1},{center + 1},{right + 1})"
                symbols = (str(atoms[left]), str(atoms[center]), str(atoms[right]))
                primitive = Primitive("angle", (left, center, right))
                specs.append(("angle", label, (left + 1, center + 1, right + 1), symbols, primitive, 180.0 / np.pi))

    for center_left, center_right in sorted(tuple(sorted(pair)) for pair in graph.bonds):
        left_neighbors = sorted(atom for atom in graph.adjacency[center_left] if atom != center_right)
        right_neighbors = sorted(atom for atom in graph.adjacency[center_right] if atom != center_left)
        for left in left_neighbors:
            for right in right_neighbors:
                if left == right:
                    continue
                label = f"D({left + 1},{center_left + 1},{center_right + 1},{right + 1})"
                symbols = (
                    str(atoms[left]),
                    str(atoms[center_left]),
                    str(atoms[center_right]),
                    str(atoms[right]),
                )
                primitive = Primitive("dihedral", (left, center_left, center_right, right))
                specs.append(
                    (
                        "dihedral",
                        label,
                        (left + 1, center_left + 1, center_right + 1, right + 1),
                        symbols,
                        primitive,
                        180.0 / np.pi,
                    )
                )

    primitives = [item[4] for item in specs]
    values = eval_primitives(primitives, coords) if primitives else np.array(())
    sigmas = _geometry_parameter_sigmas(
        primitives,
        coords,
        fit_prims=fit_prims,
        fit_u_matrix=fit_u_matrix,
        active_mask=active_mask,
        transform=transform,
        cartesian_from_parameters=cartesian_from_parameters,
        covariance=covariance,
    )
    rows: list[SemiexperimentalGeometryParameter] = []
    for idx, (kind, label, atom_indices, symbols, _primitive, angular_scale) in enumerate(specs):
        value = float(values[idx])
        sigma = sigmas[idx] if sigmas is not None and idx < len(sigmas) else None
        if kind == "bond":
            rows.append(
                SemiexperimentalGeometryParameter(
                    kind,
                    label,
                    atom_indices,
                    symbols,
                    value_angstrom=value,
                    sigma_angstrom=sigma,
                )
            )
        else:
            rows.append(
                SemiexperimentalGeometryParameter(
                    kind,
                    label,
                    atom_indices,
                    symbols,
                    value_degree=value * angular_scale,
                    sigma_degree=None if sigma is None else sigma * angular_scale,
                )
            )
    return tuple(rows)


def _geometry_parameter_sigmas(
    geometry_prims: list[Primitive],
    coords: np.ndarray,
    *,
    fit_prims: object | None,
    fit_u_matrix: np.ndarray | None,
    active_mask: np.ndarray | None,
    transform: np.ndarray | None,
    covariance: np.ndarray | None,
    cartesian_from_parameters: np.ndarray | None = None,
) -> list[float | None] | None:
    if (
        not geometry_prims
        or covariance is None
        or covariance.size == 0
    ):
        return None
    covariance = np.asarray(covariance, dtype=float)
    if cartesian_from_parameters is not None:
        dx_dr = np.asarray(cartesian_from_parameters, dtype=float)
    else:
        if (
            fit_prims is None
            or fit_u_matrix is None
            or active_mask is None
            or transform is None
            or transform.size == 0
        ):
            return None
        b_fit = np.asarray(fit_u_matrix, dtype=float).T @ b_matrix_analytic(fit_prims, coords)
        active_indices = np.where(active_mask)[0]
        dq_dr = np.zeros((b_fit.shape[0], transform.shape[1]), dtype=float)
        dq_dr[active_indices, :] = transform
        if covariance.shape != (dq_dr.shape[1], dq_dr.shape[1]):
            return None
        dx_dr = np.linalg.pinv(b_fit, rcond=1.0e-8) @ dq_dr
    if covariance.shape != (dx_dr.shape[1], dx_dr.shape[1]):
        return None
    b_geom = b_matrix_analytic(geometry_prims, coords)
    jac = b_geom @ dx_dr
    variances = np.einsum("ij,jk,ik->i", jac, covariance, jac, optimize=True)
    return [float(np.sqrt(max(value, 0.0))) if np.isfinite(value) else None for value in variances]


def kraitchman_csv_rows(rows: tuple[KraitchmanComparison, ...]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow([
        "isotopologue",
        "atom_index",
        "atom",
        "isotope_A",
        "axis",
        "kraitchman_abs_A",
        "fitted_abs_A",
        "difference_A",
    ])
    for row in rows:
        writer.writerow([
            row.isotopologue,
            row.atom_index,
            row.atom,
            row.isotope_mass_number,
            row.coordinate,
            f"{row.kraitchman_abs_angstrom:.12g}",
            f"{row.fitted_abs_angstrom:.12g}",
            f"{row.difference_angstrom:.12g}",
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


def _influence_csv(
    model: MeasurementModel | None,
    residuals: tuple[SemiexperimentalResidual, ...],
    weighted_jac: np.ndarray | None,
    weighted_residual: np.ndarray | None,
) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow([
        "row",
        "isotopologue",
        "observable",
        "residual",
        "weighted_residual",
        "chi_square_contribution",
        "leverage",
    ])
    labels = model.labels if model is not None else tuple((item.isotopologue, item.constant) for item in residuals)
    weighted = np.asarray(weighted_residual if weighted_residual is not None else (), dtype=float)
    leverage = _leverage_values(np.asarray(weighted_jac if weighted_jac is not None else np.zeros((0, 0)), dtype=float))
    for idx, item in enumerate(residuals):
        weighted_value = float(weighted[idx]) if idx < weighted.size else 0.0
        leverage_value = float(leverage[idx]) if idx < leverage.size else 0.0
        iso, obs = labels[idx] if idx < len(labels) else (item.isotopologue, item.constant)
        writer.writerow([
            idx + 1,
            iso,
            obs,
            f"{item.residual_MHz:.12g}",
            f"{weighted_value:.12g}",
            f"{weighted_value * weighted_value:.12g}",
            f"{leverage_value:.12g}",
        ])
    return stream.getvalue()


def _leverage_values(weighted_jac: np.ndarray) -> np.ndarray:
    jac = np.asarray(weighted_jac, dtype=float)
    if jac.ndim != 2 or jac.size == 0:
        return np.zeros((jac.shape[0] if jac.ndim == 2 else 0,), dtype=float)
    normal_inv = np.linalg.pinv(jac.T @ jac, rcond=1.0e-10)
    return np.einsum("ij,jk,ik->i", jac, normal_inv, jac)


def _high_correlations_csv(labels: tuple[str, ...], correlation: np.ndarray | None, threshold: float = 0.90) -> str:
    corr = np.asarray(correlation if correlation is not None else np.zeros((0, 0)), dtype=float)
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerow(["parameter_i", "parameter_j", "correlation_abs", "correlation"])
    rows = []
    for i in range(min(corr.shape[0], len(labels))):
        for j in range(i + 1, min(corr.shape[1], len(labels))):
            value = float(corr[i, j])
            abs_value = abs(value)
            if abs_value >= threshold:
                rows.append((abs_value, value, labels[i], labels[j]))
    for abs_value, value, left, right in sorted(rows, reverse=True)[:50]:
        writer.writerow([left, right, f"{abs_value:.12g}", f"{value:.12g}"])
    return stream.getvalue()


def _gic_model(
    coords: np.ndarray,
    z_numbers: np.ndarray,
    request: SemiexperimentalFitRequest | None = None,
    backend: GICForgeSEBackend | None = None,
):
    if backend is None:
        atoms = tuple(atomic_symbol(int(z)) for z in z_numbers)
        backend = _make_gicforge_backend(atoms, outdir=None)
    return backend.model(coords)


def _resolve_max_iterations(max_iter: int | None, n_optimized_parameters: int) -> int:
    if n_optimized_parameters <= 0:
        return 0
    if max_iter is not None and max_iter > 0:
        return int(max_iter)
    return max(8, 2 * int(n_optimized_parameters))


def _make_gicforge_backend(atoms: tuple[str, ...], outdir: Path | None) -> GICForgeSEBackend:
    if outdir is None:
        root = Path(tempfile.mkdtemp(prefix="merlino_se_gicforge_"))
    else:
        root = Path(outdir) / "gicforge_iterations"
        root.mkdir(parents=True, exist_ok=True)
    return GICForgeSEBackend(atoms=atoms, root=root)


def _gicforge_point_group(provout: Path) -> str:
    text = provout.read_text(encoding="utf-8", errors="replace") if provout.exists() else ""
    match = re.search(r"Point Group from symm\.f:\s*([A-Za-z0-9]+)", text)
    return match.group(1) if match else "UNKNOWN"


def _is_symmetry_refinement(previous: str, current: str) -> bool:
    previous_order = _point_group_order(previous)
    current_order = _point_group_order(current)
    return current_order > previous_order >= 1


def _point_group_order(point_group: str) -> int:
    normalized = str(point_group).strip().lower()
    explicit = {
        "c1": 1,
        "cs": 2,
        "ci": 2,
        "c2": 2,
        "c2v": 4,
        "c2h": 4,
        "d2": 4,
        "d2h": 8,
    }
    if normalized in explicit:
        return explicit[normalized]
    match = re.match(r"([cd])(\d+)([a-z]*)$", normalized)
    if not match:
        return 0
    family, order_text, suffix = match.groups()
    nfold = int(order_text)
    if family == "c":
        return 2 * nfold if suffix in {"v", "h"} else nfold
    return 4 * nfold if suffix == "h" else 2 * nfold


def _gicforge_a1_mask(labels: tuple[str, ...]) -> np.ndarray:
    irreps = []
    for label in labels:
        match = re.search(r"\birrep=([A-Za-z0-9'\"+-]+)", label)
        irreps.append(match.group(1) if match else None)
    if not any(irrep is not None for irrep in irreps):
        return np.ones(len(labels), dtype=bool)
    return np.array([irrep in {"A1", "A", "Ag", "A'"} for irrep in irreps], dtype=bool)


def _gic_model_signature(labels: tuple[str, ...]) -> tuple[int, tuple[tuple[str, int], ...], tuple[tuple[str, int], ...]]:
    irrep_counts: dict[str, int] = {}
    family_counts: dict[str, int] = {}
    for label in labels:
        irrep_match = re.search(r"\birrep=([A-Za-z0-9'\"+-]+)", label)
        irrep = irrep_match.group(1) if irrep_match else "UNK"
        family = _gic_label_family(label)
        irrep_counts[irrep] = irrep_counts.get(irrep, 0) + 1
        family_counts[family] = family_counts.get(family, 0) + 1
    return (len(labels), tuple(sorted(irrep_counts.items())), tuple(sorted(family_counts.items())))


def _validate_gic_model_signature(
    labels: tuple[str, ...],
    reference: tuple[int, tuple[tuple[str, int], ...], tuple[tuple[str, int], ...]],
) -> None:
    current = _gic_model_signature(labels)
    if current != reference:
        raise ScientificValidationError(f"GICForge coordinate model changed from {reference} to {current}")


def _gic_label_family(label: str) -> str:
    name_match = re.search(r"\bGICForge\s+([A-Za-z0-9'\"+-]+)", label)
    name = name_match.group(1) if name_match else label
    if "Str" in name:
        return "bond"
    if "Ang" in name:
        return "angle"
    if "Lin" in name:
        return "linear_bend"
    if "Tor" in name:
        return "dihedral"
    if "Oop" in name:
        return "out_of_plane"
    return "gic"


def _gic_values(prims: object, u_matrix: np.ndarray, coords: np.ndarray) -> np.ndarray:
    return u_matrix.T @ eval_primitives(prims, coords)


def _active_mask(
    labels: tuple[str, ...],
    fixed: tuple[str, ...],
    parameter_classes: tuple[ParameterClassConstraint, ...] = (),
) -> np.ndarray:
    mask = []
    fixed_l = tuple(item.lower() for item in fixed)
    fixed_classes = tuple(item for item in parameter_classes if item.mode == "fixed")
    for label in labels:
        low = label.lower()
        explicit_fixed = any(item and item in low for item in fixed_l)
        class_fixed = any(_class_matches(item, label) for item in fixed_classes)
        mask.append(not explicit_fixed and not class_fixed)
    return np.array(mask, dtype=bool)


def _gic_fixed_patterns(fixed: tuple[str, ...]) -> tuple[str, ...]:
    """Return fixed patterns that target whole GICs, not primitive coordinates."""
    return tuple(item for item in fixed if not _is_hydrogen_parameter_constraint(item) and not _primitives_from_fixed_pattern(item))


def _fixed_primitives_from_patterns(fixed: tuple[str, ...]) -> tuple[Primitive, ...]:
    primitives: list[Primitive] = []
    seen: set[tuple[str, tuple[int, ...], int]] = set()
    for item in fixed:
        for primitive in _primitives_from_fixed_pattern(item):
            key = _primitive_constraint_key(primitive)
            if key in seen:
                continue
            primitives.append(primitive)
            seen.add(key)
    return tuple(primitives)


def _hydrogen_fixed_primitives(
    atoms: list[str] | tuple[str, ...],
    available_prims: object,
    fixed: tuple[str, ...],
    *,
    coords: np.ndarray | None = None,
) -> tuple[Primitive, ...]:
    """Return a deterministic local coordinate frame for each H/D/T atom."""
    if not any(_is_hydrogen_parameter_constraint(item) for item in fixed):
        return ()
    h_atoms = {idx for idx, atom in enumerate(atoms) if str(atom).strip().upper() in {"H", "D", "T"}}
    if not h_atoms:
        return ()
    supported = {"bond", "angle", "dihedral", "out_of_plane", "linear_bend"}
    prims = [
        primitive
        for primitive in _constraint_primitive_pool(atoms, available_prims, coords)
        if primitive.kind in supported
    ]
    adjacency = _bond_adjacency(prims)
    primitives: list[Primitive] = []
    seen: set[tuple[str, tuple[int, ...], int]] = set()

    def add(primitive: Primitive | None) -> None:
        if primitive is None:
            return
        key = _primitive_constraint_key(primitive)
        if key in seen:
            return
        primitives.append(primitive)
        seen.add(key)

    for h_atom in sorted(h_atoms):
        anchors = sorted(atom for atom in adjacency.get(h_atom, ()) if atom not in h_atoms)
        if not anchors:
            # Last-resort fallback for unusual inputs: keep only directly
            # available primitives, rather than silently ignoring the H atom.
            for primitive in _hydrogen_fallback_primitives(prims, h_atom):
                add(primitive)
            continue
        anchor = anchors[0]
        add(_bond_primitive(prims, h_atom, anchor))
        linear_pair = _hydrogen_linear_pair(prims, h_atom, anchor, h_atoms)
        if linear_pair:
            for primitive in linear_pair:
                add(primitive)
            continue
        first_angle = _hydrogen_angle_primitive(prims, h_atom, anchor, h_atoms)
        add(first_angle)
        orientation = _hydrogen_orientation_primitive(prims, h_atom, anchor, h_atoms)
        if orientation is None:
            orientation = _hydrogen_angle_primitive(
                prims,
                h_atom,
                anchor,
                h_atoms,
                exclude={_primitive_constraint_key(first_angle)} if first_angle is not None else set(),
            )
        add(orientation)
    return tuple(primitives)


def _constraint_primitive_pool(
    atoms: list[str] | tuple[str, ...],
    available_prims: object,
    coords: np.ndarray | None,
) -> tuple[Primitive, ...]:
    primitives: list[Primitive] = []
    seen: set[tuple[str, tuple[int, ...], int]] = set()

    def add(primitive: Primitive) -> None:
        key = _primitive_constraint_key(primitive)
        if key in seen:
            return
        primitives.append(primitive)
        seen.add(key)

    if coords is not None:
        try:
            z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
            _continuous, graph, _ringset, _synthons, _aromaticity = build_topology_objects(
                np.asarray(coords, dtype=float),
                z_numbers,
            )
            for primitive in build_primitives(graph, np.asarray(coords, dtype=float)):
                add(primitive)
        except Exception:
            pass
    for primitive in available_prims:
        add(primitive)
    return tuple(primitives)


def _bond_adjacency(prims: list[Primitive]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {}
    for primitive in prims:
        if primitive.kind != "bond" or len(primitive.atoms) != 2:
            continue
        i, j = primitive.atoms
        adjacency.setdefault(i, set()).add(j)
        adjacency.setdefault(j, set()).add(i)
    return adjacency


def _bond_primitive(prims: list[Primitive], atom_a: int, atom_b: int) -> Primitive | None:
    wanted = {atom_a, atom_b}
    return next((primitive for primitive in prims if primitive.kind == "bond" and set(primitive.atoms) == wanted), None)


def _hydrogen_linear_pair(
    prims: list[Primitive],
    h_atom: int,
    anchor: int,
    h_atoms: set[int],
) -> tuple[Primitive, ...]:
    groups: dict[tuple[int, int, int], list[Primitive]] = {}
    for primitive in prims:
        if primitive.kind != "linear_bend" or len(primitive.atoms) != 3:
            continue
        i, j, k = primitive.atoms
        if j != anchor or h_atom not in {i, k}:
            continue
        other = k if i == h_atom else i
        key = (1 if other in h_atoms else 0, other, min(i, k))
        groups.setdefault(key, []).append(primitive)
    for _key, items in sorted(groups.items()):
        modes = {primitive.mode: primitive for primitive in items}
        if -1 in modes and -2 in modes:
            return (modes[-1], modes[-2])
    return ()


def _hydrogen_angle_primitive(
    prims: list[Primitive],
    h_atom: int,
    anchor: int,
    h_atoms: set[int],
    *,
    exclude: set[tuple[str, tuple[int, ...], int]] | None = None,
) -> Primitive | None:
    excluded = exclude or set()
    candidates = []
    for primitive in prims:
        if primitive.kind != "angle" or len(primitive.atoms) != 3:
            continue
        i, j, k = primitive.atoms
        if j != anchor or h_atom not in {i, k}:
            continue
        if _primitive_constraint_key(primitive) in excluded:
            continue
        other = k if i == h_atom else i
        candidates.append((1 if other in h_atoms else 0, other, _primitive_constraint_key(primitive), primitive))
    return min(candidates, default=(None, None, None, None))[3]


def _hydrogen_orientation_primitive(
    prims: list[Primitive],
    h_atom: int,
    anchor: int,
    h_atoms: set[int],
) -> Primitive | None:
    dihedrals = []
    for primitive in prims:
        if primitive.kind != "dihedral" or len(primitive.atoms) != 4:
            continue
        atoms = primitive.atoms
        terminal_h = (atoms[0] == h_atom and atoms[1] == anchor) or (atoms[3] == h_atom and atoms[2] == anchor)
        if not terminal_h:
            continue
        other_h_count = sum(1 for atom in atoms if atom != h_atom and atom in h_atoms)
        dihedrals.append((other_h_count, _primitive_constraint_key(primitive), primitive))
    if dihedrals:
        return min(dihedrals)[2]
    oops = []
    for primitive in prims:
        if primitive.kind != "out_of_plane" or len(primitive.atoms) != 4:
            continue
        atoms = primitive.atoms
        if atoms[0] != h_atom or atoms[1] != anchor:
            continue
        other_h_count = sum(1 for atom in atoms[2:] if atom in h_atoms)
        oops.append((other_h_count, _primitive_constraint_key(primitive), primitive))
    return min(oops, default=(None, None, None))[2]


def _hydrogen_fallback_primitives(prims: list[Primitive], h_atom: int) -> tuple[Primitive, ...]:
    candidates = [primitive for primitive in prims if h_atom in primitive.atoms]
    candidates.sort(key=lambda primitive: (
        {"bond": 0, "angle": 1, "linear_bend": 2, "dihedral": 3, "out_of_plane": 4}.get(primitive.kind, 9),
        _primitive_constraint_key(primitive),
    ))
    return tuple(candidates[:3])


def _is_hydrogen_parameter_constraint(item: str) -> bool:
    text = str(item).strip().lower().replace("-", "_").replace(" ", "_")
    return text in {
        HYDROGEN_PARAMETER_CONSTRAINT,
        "hydrogen_parameters",
        "hydrogen_primitives",
        "all_hydrogen_parameters",
        "all_hydrogen_primitives",
        "@hydrogen",
    }


def _merge_primitives(*groups: tuple[Primitive, ...]) -> tuple[Primitive, ...]:
    primitives: list[Primitive] = []
    seen: set[tuple[str, tuple[int, ...], int]] = set()
    for group in groups:
        for primitive in group:
            key = _primitive_constraint_key(primitive)
            if key in seen:
                continue
            primitives.append(primitive)
            seen.add(key)
    return tuple(primitives)


def _symmetry_expanded_fixed_primitives(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    available_prims: object,
    fixed_primitives: tuple[Primitive, ...],
) -> tuple[Primitive, ...]:
    """Expand fixed primitive constraints to the full molecular symmetry orbit."""
    if not fixed_primitives:
        return ()
    try:
        z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
        symbols = [atomic_symbol(int(z)) for z in z_numbers]
        oriented = orient_coords(coords, weights=z_numbers)
        _elements, _classes, permutations = symmetry_elements_from_geometry(
            symbols,
            oriented,
            tol=GIC_SYMM_TOL,
            max_n=6,
            tol_H=GIC_SYMM_TOL,
            ignore_isotopes=True,
            auto_max_n=True,
            inertia_tol=GIC_SYMM_INERTIA_TOL,
        )
    except Exception:
        return fixed_primitives
    if not permutations:
        return fixed_primitives

    basis = list(available_prims)
    available_by_key: dict[tuple[str, tuple[int, ...], int], Primitive] = {}
    for primitive in basis:
        available_by_key.setdefault(_primitive_constraint_key(primitive), primitive)
    basis_keys = set(available_by_key)
    for primitive in fixed_primitives:
        key = _primitive_constraint_key(primitive)
        if key not in basis_keys:
            basis_keys.add(key)
            basis.append(primitive)

    basis_positions: dict[tuple[str, tuple[int, ...], int], int] = {}
    for idx, primitive in enumerate(basis):
        basis_positions.setdefault(_primitive_constraint_key(primitive), idx)

    expanded: list[Primitive] = []
    seen: set[tuple[str, tuple[int, ...], int]] = set()

    def add(primitive: Primitive) -> None:
        key = _primitive_constraint_key(primitive)
        if key in seen:
            return
        expanded.append(available_by_key.get(key, primitive))
        seen.add(key)

    for primitive in fixed_primitives:
        add(primitive)
        seed_key = _primitive_constraint_key(primitive)
        seed_index = basis_positions.get(seed_key)
        for mapping in permutations:
            mapped = _map_primitive_by_atoms(primitive, mapping)
            mapped_key = _primitive_constraint_key(mapped)
            candidate = available_by_key.get(mapped_key, mapped)
            if seed_index is not None:
                try:
                    perm_idx, _sign = primitive_permutation(basis, mapping)
                    permuted = basis[perm_idx[seed_index]]
                    if permuted.kind == primitive.kind and _primitive_constraint_key(permuted) == mapped_key:
                        candidate = available_by_key.get(mapped_key, permuted)
                except Exception:
                    pass
            add(candidate)
    return tuple(expanded)


def _map_primitive_by_atoms(primitive: Primitive, atom_map: object) -> Primitive:
    mapped_atoms = tuple(int(atom_map[atom]) for atom in primitive.atoms)
    return Primitive(primitive.kind, mapped_atoms, mode=primitive.mode, ref=primitive.ref)


def _primitives_from_fixed_pattern(pattern: str) -> tuple[Primitive, ...]:
    text = str(pattern).strip().lower()
    match = re.match(r"^(bond|angle|dihedral|out_of_plane|linear_bend)\(([^)]*)", text)
    if not match:
        return ()
    kind, args_text = match.groups()
    args = [part.strip() for part in re.split(r"[,;]", args_text) if part.strip()]
    atoms: list[int] = []
    mode: int | None = None
    for arg in args:
        if arg.startswith("mode="):
            try:
                mode = int(arg.split("=", 1)[1])
            except ValueError:
                return ()
            continue
        try:
            atoms.append(int(arg) - 1)
        except ValueError:
            return ()
    if any(atom < 0 for atom in atoms):
        return ()
    if kind == "bond" and len(atoms) >= 2:
        return (Primitive("bond", tuple(atoms[:2])),)
    if kind == "angle" and len(atoms) >= 3:
        return (Primitive("angle", tuple(atoms[:3])),)
    if kind == "dihedral" and len(atoms) >= 4:
        return (Primitive("dihedral", tuple(atoms[:4])),)
    if kind == "out_of_plane" and len(atoms) >= 4:
        return (Primitive("out_of_plane", tuple(atoms[:4])),)
    if kind == "linear_bend" and len(atoms) >= 3:
        if mode in {-1, -2}:
            return (Primitive("linear_bend", tuple(atoms[:3]), mode=mode),)
        return (
            Primitive("linear_bend", tuple(atoms[:3]), mode=-1),
            Primitive("linear_bend", tuple(atoms[:3]), mode=-2),
        )
    return ()


def _primitive_constraint_key(primitive: Primitive) -> tuple[str, tuple[int, ...], int]:
    atoms = primitive.atoms
    if primitive.kind == "bond" and len(atoms) == 2:
        atoms = tuple(sorted(atoms))
    elif primitive.kind in {"angle", "dihedral"}:
        reverse = tuple(reversed(atoms))
        atoms = min(atoms, reverse)
    elif primitive.kind == "out_of_plane" and len(atoms) == 4:
        atoms = (atoms[1], *tuple(sorted((atoms[0], atoms[2], atoms[3]))))
    elif primitive.kind == "linear_bend" and len(atoms) == 3:
        atoms = (atoms[1], *tuple(sorted((atoms[0], atoms[2]))))
    return (primitive.kind, tuple(atoms), primitive.mode)


def _primitive_constrained_transform(
    coords: np.ndarray,
    prims: object,
    u_matrix: np.ndarray,
    active_mask: np.ndarray,
    transform: np.ndarray,
    names: tuple[str, ...],
    fixed_primitives: tuple[Primitive, ...],
    *,
    cartesian_from_q: np.ndarray | None = None,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Project reduced GIC increments onto the null space of fixed primitives."""
    if not fixed_primitives or transform.size == 0:
        return transform, names
    active_indices = np.where(active_mask)[0]
    if not len(active_indices):
        return transform, names
    if cartesian_from_q is None:
        cartesian_from_q = _gic_cartesian_projector(prims, u_matrix, coords)
    if cartesian_from_q.size == 0:
        return transform, names
    b_fixed = b_matrix_analytic(list(fixed_primitives), coords)
    constraints_active = (b_fixed @ cartesian_from_q)[:, active_indices]
    constraints_reduced = constraints_active @ transform
    constraints_reduced = _independent_rows_incremental(constraints_reduced)
    null = _nullspace(constraints_reduced)
    constrained = transform @ null
    if constrained.shape[1] == len(names):
        return constrained, names
    return constrained, tuple(f"constrained_{idx:03d}" for idx in range(1, constrained.shape[1] + 1))


def _primitive_constrained_cartesian_transform(
    coords: np.ndarray,
    cartesian_from_q: np.ndarray,
    active_mask: np.ndarray,
    transform: np.ndarray,
    names: tuple[str, ...],
    fixed_primitives: tuple[Primitive, ...],
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Project reduced Cartesian-basis increments onto fixed primitive constraints."""
    if not fixed_primitives or transform.size == 0:
        return transform, names
    active_indices = np.where(active_mask)[0]
    if not len(active_indices):
        return transform, names
    basis = np.asarray(cartesian_from_q, dtype=float)
    if basis.size == 0:
        return transform, names
    b_fixed = b_matrix_analytic(list(fixed_primitives), coords)
    constraints_active = (b_fixed @ basis)[:, active_indices]
    constraints_reduced = constraints_active @ transform
    constraints_reduced = _independent_rows_incremental(constraints_reduced)
    null = _nullspace(constraints_reduced)
    constrained = transform @ null
    if constrained.shape[1] == len(names):
        return constrained, names
    return constrained, tuple(f"constrained_{idx:03d}" for idx in range(1, constrained.shape[1] + 1))


def _cartesian_from_reduced_coordinates(
    cartesian_from_q: np.ndarray,
    active_mask: np.ndarray,
    transform: np.ndarray,
) -> np.ndarray:
    active_indices = np.where(active_mask)[0]
    if transform.size == 0 or not len(active_indices):
        return np.zeros((np.asarray(cartesian_from_q).shape[0], 0), dtype=float)
    return np.asarray(cartesian_from_q, dtype=float)[:, active_indices] @ transform


def _active_coordinate_jacobian(jacobian: np.ndarray, active_mask: np.ndarray) -> np.ndarray:
    return np.asarray(jacobian, dtype=float)[:, np.where(active_mask)[0]]


def _nullspace(matrix: np.ndarray) -> np.ndarray:
    mat = np.asarray(matrix, dtype=float)
    if mat.ndim != 2:
        raise ValueError("Null-space input must be a matrix")
    ncols = mat.shape[1]
    if ncols == 0:
        return np.zeros((0, 0), dtype=float)
    if mat.size == 0:
        return np.eye(ncols, dtype=float)
    _u, singular, vh = np.linalg.svd(mat, full_matrices=True)
    if not singular.size:
        return np.eye(ncols, dtype=float)
    tol = max(mat.shape) * np.finfo(float).eps * float(singular[0])
    rank = int(np.sum(singular > tol))
    return vh[rank:, :].T.copy()


def _independent_rows_incremental(matrix: np.ndarray) -> np.ndarray:
    mat = np.asarray(matrix, dtype=float)
    if mat.ndim != 2 or mat.size == 0:
        return mat
    row_norms = np.linalg.norm(mat, axis=1)
    scale = float(np.max(row_norms)) if row_norms.size else 0.0
    tol = max(mat.shape) * np.finfo(float).eps * max(scale, 1.0) * 100.0
    basis: list[np.ndarray] = []
    selected: list[int] = []
    for idx, row in enumerate(mat):
        residual = row.astype(float, copy=True)
        for vector in basis:
            residual -= vector * float(vector @ residual)
        norm = float(np.linalg.norm(residual))
        if norm > tol:
            basis.append(residual / norm)
            selected.append(idx)
    return mat[selected, :] if selected else np.zeros((0, mat.shape[1]), dtype=float)


def _incremental_column_rank(matrix: np.ndarray) -> int:
    mat = np.asarray(matrix, dtype=float)
    if mat.ndim != 2 or mat.size == 0:
        return 0
    col_norms = np.linalg.norm(mat, axis=0)
    scale = float(np.max(col_norms)) if col_norms.size else 0.0
    tol = max(mat.shape) * np.finfo(float).eps * max(scale, 1.0) * 100.0
    basis: list[np.ndarray] = []
    for col in range(mat.shape[1]):
        residual = mat[:, col].astype(float, copy=True)
        for vector in basis:
            residual -= vector * float(vector @ residual)
        norm = float(np.linalg.norm(residual))
        if norm > tol:
            basis.append(residual / norm)
    return len(basis)


def _auto_pruned_active_mask(labels: tuple[str, ...], patterns: tuple[str, ...]) -> np.ndarray:
    if not patterns:
        return np.ones(len(labels), dtype=bool)
    lowered = tuple(pattern.lower() for pattern in patterns)
    return np.array([not any(pattern in label.lower() for pattern in lowered) for label in labels], dtype=bool)


def _mark_auto_pruned_classes(
    labels: tuple[str, ...],
    class_by_gic: tuple[str, ...],
    patterns: tuple[str, ...],
) -> tuple[str, ...]:
    if not patterns:
        return class_by_gic
    lowered = tuple(pattern.lower() for pattern in patterns)
    classes = list(class_by_gic)
    if len(classes) < len(labels):
        classes.extend("" for _ in range(len(labels) - len(classes)))
    for idx, label in enumerate(labels):
        if any(pattern in label.lower() for pattern in lowered):
            classes[idx] = "auto_pruned_weak"
    return tuple(classes)


def _weak_parameter_patterns(
    names: tuple[str, ...],
    weighted_jac: np.ndarray,
    condition_target: float,
) -> tuple[str, ...]:
    if weighted_jac.size == 0 or weighted_jac.shape[1] <= 1 or condition_target <= 0.0:
        return ()
    remaining = list(range(weighted_jac.shape[1]))
    pruned: list[str] = []
    while len(remaining) > 1:
        current = weighted_jac[:, remaining]
        conditioning = rank_condition(current)
        if np.isfinite(conditioning.condition_number) and conditioning.condition_number <= condition_target:
            break
        best: tuple[float, int] | None = None
        for col in remaining:
            trial = [item for item in remaining if item != col]
            trial_condition = rank_condition(weighted_jac[:, trial]).condition_number
            if not np.isfinite(trial_condition):
                continue
            score = (trial_condition, col)
            if best is None or score < best:
                best = score
        if best is None or best[0] >= conditioning.condition_number:
            break
        removed = best[1]
        pattern = _parameter_prune_pattern(names[removed])
        if pattern:
            pruned.append(pattern)
        remaining.remove(removed)
    return tuple(pruned)


def _parameter_prune_pattern(name: str) -> str:
    parts = str(name).split()
    if len(parts) >= 3 and re.match(r"^[A-Z][0-9][A-Za-z]+[0-9]+$", parts[2]):
        return parts[2]
    if parts:
        return parts[0]
    return str(name)


def _parameter_class_transform(
    labels: tuple[str, ...],
    active_mask: np.ndarray,
    parameter_classes: tuple[ParameterClassConstraint, ...],
) -> tuple[np.ndarray, tuple[str, ...], tuple[str, ...]]:
    active_indices = np.where(active_mask)[0]
    class_by_gic = [
        next((item.name for item in parameter_classes if item.mode == "fixed" and _class_matches(item, label)), "")
        for label in labels
    ]
    if not len(active_indices):
        return np.zeros((0, 0), dtype=float), (), tuple(class_by_gic)
    columns: list[np.ndarray] = []
    names: list[str] = []
    shared_classes = tuple(item for item in parameter_classes if item.mode == "shared")
    assigned = np.zeros(len(active_indices), dtype=bool)
    for parameter_class in shared_classes:
        local = [pos for pos, idx in enumerate(active_indices) if _class_matches(parameter_class, labels[idx])]
        if not local:
            continue
        col = np.zeros(len(active_indices), dtype=float)
        for pos in local:
            col[pos] = 1.0
            assigned[pos] = True
            class_by_gic[active_indices[pos]] = parameter_class.name
        columns.append(col)
        names.append(parameter_class.name)
    for pos, idx in enumerate(active_indices):
        if assigned[pos]:
            continue
        col = np.zeros(len(active_indices), dtype=float)
        col[pos] = 1.0
        columns.append(col)
        names.append(labels[idx])
    return np.column_stack(columns), tuple(names), tuple(class_by_gic)


def _reduced_parameter_scales(labels: tuple[str, ...], active_mask: np.ndarray, transform: np.ndarray) -> np.ndarray:
    if transform.size == 0:
        return np.ones(transform.shape[1], dtype=float)
    active_indices = np.where(active_mask)[0]
    gic_scales = np.array([_gic_coordinate_scale(labels[idx]) for idx in active_indices], dtype=float)
    scales = np.ones(transform.shape[1], dtype=float)
    for col in range(transform.shape[1]):
        weights = np.abs(transform[:, col])
        weight_sum = float(np.sum(weights))
        if weight_sum > 0.0:
            scales[col] = float(np.sum(weights * gic_scales) / weight_sum)
    return np.clip(scales, 0.05, 2.0)


def _gic_coordinate_scale(label: str) -> float:
    low = label.lower()
    if "bond(" in low or "str" in low:
        return 1.0
    if any(token in low for token in ("angle(", "dihedral(", "out_of_plane(", "linear_bend(", "ang", "dih", "oop", "lin", "pck")):
        return 0.5
    return 1.0


def _class_matches(parameter_class: ParameterClassConstraint, label: str) -> bool:
    low = label.lower()
    return any(pattern.lower() in low for pattern in parameter_class.patterns)


def _jacobian_constants_wrt_gics(
    atoms: list[str],
    coords: np.ndarray,
    request: SemiexperimentalFitRequest,
    prims: object,
    u_matrix: np.ndarray,
    active_mask: np.ndarray,
    labels: tuple[str, ...],
    measurement_model: "MeasurementModel",
    *,
    step: float,
    cartesian_from_q: np.ndarray | None = None,
) -> np.ndarray:
    active_indices = np.where(active_mask)[0]
    base_q = _gic_values(prims, u_matrix, coords)
    jac = np.zeros((len(measurement_model.observed), len(active_indices)), dtype=float)
    if cartesian_from_q is None:
        cartesian_from_q = _gic_cartesian_projector(prims, u_matrix, coords)
    analytic = _analytic_measurement_jacobian_wrt_gics(
        atoms,
        coords,
        request,
        labels,
        measurement_model,
        cartesian_from_q,
    )
    if analytic is not None and analytic.shape == (len(measurement_model.observed), len(base_q)):
        return analytic[:, active_indices]
    return _finite_difference_measurement_jacobian_wrt_gics(
        atoms,
        coords,
        request,
        prims,
        u_matrix,
        active_indices,
        labels,
        measurement_model,
        step=step,
        cartesian_from_q=cartesian_from_q,
    )


def _analytic_measurement_jacobian_wrt_gics(
    atoms: list[str],
    coords: np.ndarray,
    request: SemiexperimentalFitRequest,
    labels: tuple[str, ...],
    measurement_model: "MeasurementModel",
    cartesian_from_q: np.ndarray,
) -> np.ndarray | None:
    if measurement_model.observable == "moments":
        cartesian = _moments_cartesian_jacobian(atoms, coords, request.observations)
        selected = _select_raw_components(cartesian, MOMENT_COMPONENTS, measurement_model.components)
    elif measurement_model.observable == "rotational_constants":
        cartesian = _rotational_constants_cartesian_jacobian(atoms, coords, request.observations)
        selected = _select_raw_components(cartesian, ROTATIONAL_COMPONENTS, measurement_model.components)
    else:
        return None
    gic_jac = selected @ cartesian_from_q
    predicate = _predicate_jacobian(request.qm_predicates, labels, cartesian_from_q.shape[1])
    if predicate.size:
        return np.vstack([gic_jac, predicate])
    return gic_jac


def _finite_difference_measurement_jacobian_wrt_gics(
    atoms: list[str],
    coords: np.ndarray,
    request: SemiexperimentalFitRequest,
    prims: object,
    u_matrix: np.ndarray,
    active_indices: np.ndarray,
    labels: tuple[str, ...],
    measurement_model: "MeasurementModel",
    *,
    step: float,
    cartesian_from_q: np.ndarray,
) -> np.ndarray:
    base_q = _gic_values(prims, u_matrix, coords)
    jac = np.zeros((len(measurement_model.observed), len(active_indices)), dtype=float)

    def column(idx: int) -> np.ndarray:
        dq = np.zeros_like(base_q)
        dq[idx] = step
        plus = _displace_along_gics(coords, prims, u_matrix, dq, cartesian_from_q=cartesian_from_q)
        plus_q = _gic_values(prims, u_matrix, plus)
        dq[idx] = -step
        minus = _displace_along_gics(coords, prims, u_matrix, dq, cartesian_from_q=cartesian_from_q)
        minus_q = _gic_values(prims, u_matrix, minus)
        return (
            _measurement_vector(atoms, plus, request, plus_q, labels, measurement_model)
            - _measurement_vector(atoms, minus, request, minus_q, labels, measurement_model)
        ) / (2.0 * step)

    max_workers = min(len(active_indices), max(1, (os.cpu_count() or 1)))
    if max_workers <= 1 or len(active_indices) < 4:
        for col, idx in enumerate(active_indices):
            jac[:, col] = column(int(idx))
        return jac
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for col, values in enumerate(executor.map(column, [int(idx) for idx in active_indices])):
            jac[:, col] = values
    return jac


def _jacobian_constants_wrt_cartesian_basis(
    atoms: list[str],
    coords: np.ndarray,
    request: SemiexperimentalFitRequest,
    labels: tuple[str, ...],
    measurement_model: "MeasurementModel",
    cartesian_from_q: np.ndarray,
) -> np.ndarray:
    if measurement_model.observable == "moments":
        cartesian = _moments_cartesian_jacobian(atoms, coords, request.observations)
        selected = _select_raw_components(cartesian, MOMENT_COMPONENTS, measurement_model.components)
    elif measurement_model.observable == "rotational_constants":
        cartesian = _rotational_constants_cartesian_jacobian(atoms, coords, request.observations)
        selected = _select_raw_components(cartesian, ROTATIONAL_COMPONENTS, measurement_model.components)
    else:
        raise ScientificValidationError(f"Unsupported observable for Cartesian-basis SEfit: {measurement_model.observable}")
    jac = selected @ cartesian_from_q
    predicate = _predicate_jacobian(request.qm_predicates, labels, cartesian_from_q.shape[1])
    if predicate.size:
        return np.vstack([jac, predicate])
    return jac


def _moments_cartesian_jacobian(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
) -> np.ndarray:
    rows = []
    for obs in observations:
        isotopes = _isotopes_for_observation(atoms, obs)
        _moments, jac = _principal_moments_and_cartesian_jacobian(atoms, coords, isotopes)
        rows.append(jac)
    return np.vstack(rows) if rows else np.zeros((0, np.asarray(coords).size), dtype=float)


def _rotational_constants_cartesian_jacobian(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
) -> np.ndarray:
    rows = []
    for obs in observations:
        isotopes = _isotopes_for_observation(atoms, obs)
        moments, moment_jac = _principal_moments_and_cartesian_jacobian(atoms, coords, isotopes)
        factors = np.zeros(3, dtype=float)
        positive = moments > 0.0
        factors[positive] = -ROTCONST_TO_MOMENT / (moments[positive] * moments[positive])
        rows.append(factors[:, None] * moment_jac)
    return np.vstack(rows) if rows else np.zeros((0, np.asarray(coords).size), dtype=float)


def _principal_moments_and_cartesian_jacobian(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    isotopes: list[int | None],
) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(coords, dtype=float)
    structure = Structure.from_atoms_coords(list(atoms), [tuple(row) for row in arr], isotopes=isotopes)
    masses = np.asarray(structure.mass_isotope, dtype=float)
    total_mass = float(np.sum(masses))
    centered = arr - (masses[:, None] * arr).sum(axis=0) / total_mass
    inertia = np.zeros((3, 3), dtype=float)
    eye = np.eye(3)
    for mass, xyz in zip(masses, centered):
        inertia += mass * ((xyz @ xyz) * eye - np.outer(xyz, xyz))
    moments, axes = np.linalg.eigh(inertia)
    jac = np.zeros((3, arr.size), dtype=float)
    for atom_idx, (mass, xyz) in enumerate(zip(masses, centered)):
        for axis_idx in range(3):
            unit = eye[axis_idx]
            derivative = mass * (2.0 * xyz[axis_idx] * eye - np.outer(unit, xyz) - np.outer(xyz, unit))
            col = 3 * atom_idx + axis_idx
            for moment_idx in range(3):
                vector = axes[:, moment_idx]
                jac[moment_idx, col] = float(vector @ derivative @ vector)
    return moments, jac


def _predicate_jacobian(
    predicates: tuple[QMParameterPredicate, ...],
    labels: tuple[str, ...],
    n_q: int,
) -> np.ndarray:
    rows = []
    for predicate in predicates:
        for idx in _predicate_indices(predicate, labels):
            row = np.zeros(n_q, dtype=float)
            row[idx] = 1.0
            rows.append(row)
    return np.vstack(rows) if rows else np.zeros((0, n_q), dtype=float)


def _gic_cartesian_projector(prims: object, u_matrix: np.ndarray, coords: np.ndarray) -> np.ndarray:
    bq = u_matrix.T @ b_matrix_analytic(prims, coords)
    return np.linalg.pinv(bq, rcond=1.0e-8)


def _gic_projector_state(
    prims: object,
    u_matrix: np.ndarray,
    coords: np.ndarray,
    q_values: np.ndarray,
) -> GICProjectorState:
    return GICProjectorState(
        coords=np.asarray(coords, dtype=float).copy(),
        q_values=np.asarray(q_values, dtype=float).copy(),
        cartesian_from_q=_gic_cartesian_projector(prims, u_matrix, coords),
    )


def _secant_projector_update(
    cartesian_from_q: np.ndarray,
    previous_coords: np.ndarray,
    previous_q: np.ndarray,
    current_coords: np.ndarray,
    current_q: np.ndarray,
) -> SecantProjectorUpdate:
    q_delta = np.asarray(current_q, dtype=float) - np.asarray(previous_q, dtype=float)
    x_delta = np.asarray(current_coords, dtype=float).reshape(-1) - np.asarray(previous_coords, dtype=float).reshape(-1)
    if cartesian_from_q.shape != (x_delta.size, q_delta.size):
        return SecantProjectorUpdate(None, float("inf"), False)
    denom = float(q_delta @ q_delta)
    if denom <= 1.0e-24 or not np.isfinite(denom):
        return SecantProjectorUpdate(None, float("inf"), False)
    predicted_x_delta = cartesian_from_q @ q_delta
    residual = x_delta - predicted_x_delta
    x_norm = float(np.linalg.norm(x_delta))
    relative_error = float(np.linalg.norm(residual) / max(x_norm, 1.0e-12))
    update = np.outer(residual, q_delta) / denom
    updated = cartesian_from_q + update
    if not np.all(np.isfinite(updated)):
        return SecantProjectorUpdate(None, relative_error, False)
    update_norm = float(np.linalg.norm(update))
    projector_norm = float(np.linalg.norm(cartesian_from_q))
    if relative_error > 0.75 or update_norm > max(0.75 * projector_norm, 1.0e-8):
        return SecantProjectorUpdate(None, relative_error, False)
    return SecantProjectorUpdate(updated, relative_error, True)


def _displace_along_gics(
    coords: np.ndarray,
    prims: object,
    u_matrix: np.ndarray,
    dq: np.ndarray,
    *,
    cartesian_from_q: np.ndarray | None = None,
) -> np.ndarray:
    if cartesian_from_q is None:
        cartesian_from_q = _gic_cartesian_projector(prims, u_matrix, coords)
    dx = cartesian_from_q @ dq
    return coords + dx.reshape(coords.shape)


def _displace_along_cartesian_basis(
    coords: np.ndarray,
    cartesian_from_q: np.ndarray,
    dq: np.ndarray,
) -> np.ndarray:
    dx = np.asarray(cartesian_from_q, dtype=float) @ np.asarray(dq, dtype=float)
    return np.asarray(coords, dtype=float) + dx.reshape(np.asarray(coords).shape)


def _line_search_update(
    atoms: list[str],
    coords: np.ndarray,
    request: SemiexperimentalFitRequest,
    labels: tuple[str, ...],
    measurement_model: "MeasurementModel",
    prims: object,
    u_matrix: np.ndarray,
    dq: np.ndarray,
    *,
    current_objective: float,
    base_q: np.ndarray,
    cartesian_from_q: np.ndarray | None = None,
    weighted_residual: np.ndarray,
    jac_weighted: np.ndarray,
    reduced_step: np.ndarray,
) -> LineSearchResult:
    observed = measurement_model.observed
    sqrt_weights = np.sqrt(measurement_model.weights)
    best = LineSearchResult(
        coords=coords,
        q_values=base_q,
        objective=current_objective,
        accepted=False,
        actual_reduction=0.0,
        predicted_reduction=0.0,
        ratio=0.0,
        scale=0.0,
    )
    if cartesian_from_q is None:
        cartesian_from_q = _gic_cartesian_projector(prims, u_matrix, coords)
    for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625):
        predicted_reduction = _predicted_reduction(
            weighted_residual,
            jac_weighted,
            reduced_step,
            scale=scale,
            current_objective=current_objective,
        )
        candidate = _displace_along_gics(coords, prims, u_matrix, scale * dq, cartesian_from_q=cartesian_from_q)
        try:
            q_candidate = _gic_values(prims, u_matrix, candidate)
            calc = _measurement_vector(atoms, candidate, request, q_candidate, labels, measurement_model)
        except Exception:
            continue
        if calc.shape != observed.shape:
            continue
        residual = (observed - calc) * sqrt_weights
        candidate_objective = objective(residual)
        if not np.isfinite(candidate_objective):
            continue
        actual_reduction = current_objective - candidate_objective
        ratio = actual_reduction / predicted_reduction if predicted_reduction > 0.0 else 0.0
        if actual_reduction > 0.0 and candidate_objective < best.objective:
            best = LineSearchResult(
                coords=candidate,
                q_values=q_candidate,
                objective=candidate_objective,
                accepted=True,
                actual_reduction=float(actual_reduction),
                predicted_reduction=float(max(predicted_reduction, 0.0)),
                ratio=float(ratio),
                scale=float(scale),
            )
            if predicted_reduction <= 0.0 or actual_reduction >= 1.0e-4 * predicted_reduction:
                break
    return best


def _line_search_update_cartesian_basis(
    atoms: list[str],
    coords: np.ndarray,
    request: SemiexperimentalFitRequest,
    labels: tuple[str, ...],
    measurement_model: "MeasurementModel",
    mode_model: CartesianCoordinateModel,
    dq: np.ndarray,
    *,
    current_objective: float,
    base_q: np.ndarray,
    weighted_residual: np.ndarray,
    jac_weighted: np.ndarray,
    reduced_step: np.ndarray,
) -> LineSearchResult:
    observed = measurement_model.observed
    sqrt_weights = np.sqrt(measurement_model.weights)
    best = LineSearchResult(
        coords=coords,
        q_values=base_q,
        objective=current_objective,
        accepted=False,
        actual_reduction=0.0,
        predicted_reduction=0.0,
        ratio=0.0,
        scale=0.0,
    )
    for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625):
        predicted_reduction = _predicted_reduction(
            weighted_residual,
            jac_weighted,
            reduced_step,
            scale=scale,
            current_objective=current_objective,
        )
        candidate = _displace_along_cartesian_basis(coords, mode_model.cartesian_from_q, scale * dq)
        try:
            q_candidate = mode_model.values(candidate)
            calc = _measurement_vector(atoms, candidate, request, q_candidate, labels, measurement_model)
        except Exception:
            continue
        if calc.shape != observed.shape:
            continue
        residual = (observed - calc) * sqrt_weights
        candidate_objective = objective(residual)
        if not np.isfinite(candidate_objective):
            continue
        actual_reduction = current_objective - candidate_objective
        ratio = actual_reduction / predicted_reduction if predicted_reduction > 0.0 else 0.0
        if actual_reduction > 0.0 and candidate_objective < best.objective:
            best = LineSearchResult(
                coords=candidate,
                q_values=q_candidate,
                objective=candidate_objective,
                accepted=True,
                actual_reduction=float(actual_reduction),
                predicted_reduction=float(max(predicted_reduction, 0.0)),
                ratio=float(ratio),
                scale=float(scale),
            )
            if predicted_reduction <= 0.0 or actual_reduction >= 1.0e-4 * predicted_reduction:
                break
    return best


def _adaptive_lm_step(
    jac_weighted: np.ndarray,
    weighted_residual: np.ndarray,
    damping: float,
    trust_radius: float,
) -> np.ndarray:
    if jac_weighted.size == 0 or jac_weighted.shape[1] == 0:
        return np.zeros((0,), dtype=float)
    try:
        step = damped_normal_step(jac_weighted, weighted_residual, damping)
    except np.linalg.LinAlgError:
        lhs = jac_weighted.T @ jac_weighted + float(damping) * np.eye(jac_weighted.shape[1])
        rhs = jac_weighted.T @ weighted_residual
        step = np.linalg.lstsq(lhs, rhs, rcond=1.0e-10)[0]
    if not np.all(np.isfinite(step)):
        step = _cauchy_step(jac_weighted, weighted_residual)
    step = limit_step(step, trust_radius)
    predicted = _predicted_reduction(
        weighted_residual,
        jac_weighted,
        step,
        scale=1.0,
        current_objective=objective(weighted_residual),
    )
    if predicted <= 0.0 or not np.all(np.isfinite(step)):
        step = limit_step(_cauchy_step(jac_weighted, weighted_residual), trust_radius)
    return step


def _cauchy_step(jac_weighted: np.ndarray, weighted_residual: np.ndarray) -> np.ndarray:
    gradient = jac_weighted.T @ weighted_residual
    if gradient.size == 0:
        return gradient
    jg = jac_weighted @ gradient
    denom = float(jg @ jg)
    if denom <= 0.0 or not np.isfinite(denom):
        norm = float(np.linalg.norm(gradient))
        return gradient / norm if norm > 0.0 else gradient
    alpha = float((gradient @ gradient) / denom)
    return alpha * gradient


def _predicted_reduction(
    weighted_residual: np.ndarray,
    jac_weighted: np.ndarray,
    reduced_step: np.ndarray,
    *,
    scale: float,
    current_objective: float,
) -> float:
    if reduced_step.size == 0:
        return 0.0
    predicted_residual = weighted_residual - float(scale) * (jac_weighted @ reduced_step)
    predicted_objective = objective(predicted_residual)
    reduction = float(current_objective - predicted_objective)
    return reduction if np.isfinite(reduction) else 0.0


def _accepted_trust_update(
    damping: float,
    trust_radius: float,
    ratio: float,
    scale: float,
    step_norm: float,
    max_step: float,
) -> tuple[float, float]:
    ratio = float(ratio) if np.isfinite(ratio) else 0.0
    if ratio < 0.25:
        new_damping = min(max(float(damping) * 4.0, 1.0e-12), 1.0e12)
        new_radius = _scaled_trust_radius(trust_radius, max_step, 0.5)
    elif ratio > 0.75 and scale >= 0.9:
        new_damping = max(float(damping) / 2.5, 1.0e-14)
        new_radius = _expanded_trust_radius(trust_radius, max_step, step_norm)
    else:
        new_damping = max(float(damping) / 1.4, 1.0e-14)
        new_radius = trust_radius
    return new_damping, new_radius


def _rejected_trust_update(damping: float, trust_radius: float, max_step: float) -> tuple[float, float]:
    return min(max(float(damping) * 6.0, 1.0e-12), 1.0e12), _scaled_trust_radius(
        trust_radius, max_step, 0.5
    )


def _scaled_trust_radius(trust_radius: float, max_step: float, scale: float) -> float:
    if max_step <= 0.0:
        return trust_radius
    current = trust_radius if trust_radius > 0.0 else max_step
    return max(float(current) * float(scale), 1.0e-7)


def _expanded_trust_radius(trust_radius: float, max_step: float, step_norm: float) -> float:
    if max_step <= 0.0:
        return trust_radius
    current = trust_radius if trust_radius > 0.0 else max_step
    proposed = max(current * 1.6, step_norm * 1.25, 1.0e-7)
    return min(float(max_step), float(proposed))


def _should_refresh_gic_model(
    line_search: LineSearchResult,
    model_age: int,
    *,
    secant_relative_error: float,
    tolerance_MHz: float,
    n_observations: int,
) -> bool:
    if model_age <= 0:
        return False
    if model_age >= 3:
        return True
    if line_search.scale < 0.5:
        return True
    if not np.isfinite(secant_relative_error) or secant_relative_error > 0.5:
        return True
    if line_search.ratio < 0.25 or line_search.ratio > 2.5:
        return True
    objective_scale = max(int(n_observations), 1) * tolerance_MHz * tolerance_MHz
    return line_search.objective <= max(10.0 * objective_scale, 1.0e-24)


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


def _moments_vector(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
) -> np.ndarray:
    values: list[float] = []
    for obs in observations:
        isotopes = _isotopes_for_observation(atoms, obs)
        structure = Structure.from_atoms_coords(list(atoms), [tuple(row) for row in coords], isotopes=isotopes)
        values.extend(principal_moments(structure, isotopic=True))
    return np.array(values, dtype=float)


def _build_measurement_model(
    request: SemiexperimentalFitRequest,
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    prims: object,
    u_matrix: np.ndarray,
    labels: tuple[str, ...],
) -> MeasurementModel:
    observable = "moments" if request.observable == "auto" else request.observable
    planar = _is_planar(coords)
    components = _select_components(request, observable, atoms, coords, prims, u_matrix, planar)
    observed = _experimental_observed_vector(request, observable, components)
    weights = _experimental_weights_vector(request, observable, components)
    row_labels: list[tuple[str, str]] = []
    for obs in request.observations:
        row_labels.extend((obs.label, comp) for comp in components)
    predicate_values, predicate_weights, predicate_labels = _predicate_observations(request.qm_predicates, labels)
    if predicate_values.size:
        observed = np.concatenate([observed, predicate_values])
        weights = np.concatenate([weights, predicate_weights])
        row_labels.extend(predicate_labels)
    return MeasurementModel(
        observable=observable,
        components=components,
        labels=tuple(row_labels),
        observed=observed,
        weights=weights,
        planar=planar,
    )


def _build_measurement_model_cartesian_basis(
    request: SemiexperimentalFitRequest,
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    labels: tuple[str, ...],
    cartesian_from_q: np.ndarray,
) -> MeasurementModel:
    observable = "moments" if request.observable == "auto" else request.observable
    planar = _is_planar(coords)
    if observable == "moments":
        components = MOMENT_COMPONENTS
    elif request.rotational_components != "auto":
        components = tuple(request.rotational_components)
    elif planar:
        components = _best_planar_rotational_pair_from_cartesian_basis(
            atoms,
            coords,
            request.observations,
            cartesian_from_q,
        )
    else:
        components = ROTATIONAL_COMPONENTS
    observed = _experimental_observed_vector(request, observable, components)
    weights = _experimental_weights_vector(request, observable, components)
    row_labels: list[tuple[str, str]] = []
    for obs in request.observations:
        row_labels.extend((obs.label, comp) for comp in components)
    predicate_values, predicate_weights, predicate_labels = _predicate_observations(request.qm_predicates, labels)
    if predicate_values.size:
        observed = np.concatenate([observed, predicate_values])
        weights = np.concatenate([weights, predicate_weights])
        row_labels.extend(predicate_labels)
    return MeasurementModel(
        observable=observable,
        components=components,
        labels=tuple(row_labels),
        observed=observed,
        weights=weights,
        planar=planar,
    )


def _measurement_vector(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    request: SemiexperimentalFitRequest,
    q_values: np.ndarray,
    labels: tuple[str, ...],
    model: MeasurementModel,
) -> np.ndarray:
    if model.observable == "moments":
        raw = _moments_vector(atoms, coords, request.observations)
        selected = _select_raw_components(raw, MOMENT_COMPONENTS, model.components)
    else:
        raw = _constants_vector(atoms, coords, request.observations)
        selected = _select_raw_components(raw, ROTATIONAL_COMPONENTS, model.components)
    predicate_values = _predicate_values(request.qm_predicates, labels, q_values)
    if predicate_values.size:
        return np.concatenate([selected, predicate_values])
    return selected


def _experimental_observed_vector(
    request: SemiexperimentalFitRequest,
    observable: str,
    components: tuple[str, ...],
) -> np.ndarray:
    if observable == "moments":
        raw = []
        for obs in request.observations:
            raw.extend(_constants_to_moments(obs.corrected.as_tuple()))
        return _select_raw_components(np.array(raw, dtype=float), MOMENT_COMPONENTS, components)
    raw = _observed_vector(request.observations)
    return _select_raw_components(raw, ROTATIONAL_COMPONENTS, components)


def _experimental_weights_vector(
    request: SemiexperimentalFitRequest,
    observable: str,
    components: tuple[str, ...],
) -> np.ndarray:
    values: list[float] = []
    for obs in request.observations:
        if observable == "moments":
            values.extend(_moment_weights(obs))
        else:
            values.extend(obs.weights.as_tuple() if obs.weights is not None else (1.0, 1.0, 1.0))
    component_names = MOMENT_COMPONENTS if observable == "moments" else ROTATIONAL_COMPONENTS
    return _select_raw_components(np.array(values, dtype=float), component_names, components)


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


def _select_components(
    request: SemiexperimentalFitRequest,
    observable: str,
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    prims: object,
    u_matrix: np.ndarray,
    planar: bool,
) -> tuple[str, ...]:
    if observable == "moments":
        return MOMENT_COMPONENTS
    if request.rotational_components != "auto":
        return tuple(request.rotational_components)
    if not planar:
        return ROTATIONAL_COMPONENTS
    return _best_planar_rotational_pair(atoms, coords, request.observations, prims, u_matrix)


def _best_planar_rotational_pair(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
    prims: object,
    u_matrix: np.ndarray,
) -> tuple[str, ...]:
    candidates = (("A", "B"), ("A", "C"), ("B", "C"))
    cartesian_from_q = _gic_cartesian_projector(prims, u_matrix, coords)
    full = _rotational_constants_cartesian_jacobian(atoms, coords, observations) @ cartesian_from_q
    best = candidates[0]
    best_score = (-1, float("inf"))
    for pair in candidates:
        subset = _select_raw_components(full, ROTATIONAL_COMPONENTS, pair)
        singular = np.linalg.svd(subset, compute_uv=False)
        rank = int(np.sum(singular > max(subset.shape) * np.finfo(float).eps * (singular[0] if singular.size else 0.0)))
        cond = float(singular[0] / singular[-1]) if singular.size and singular[-1] > 0.0 else float("inf")
        score = (rank, cond)
        if score[0] > best_score[0] or (score[0] == best_score[0] and score[1] < best_score[1]):
            best = pair
            best_score = score
    return best


def _best_planar_rotational_pair_from_cartesian_basis(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
    cartesian_from_q: np.ndarray,
) -> tuple[str, ...]:
    candidates = (("A", "B"), ("A", "C"), ("B", "C"))
    full = _rotational_constants_cartesian_jacobian(atoms, coords, observations) @ cartesian_from_q
    best = candidates[0]
    best_score = (-1, float("inf"))
    for pair in candidates:
        subset = _select_raw_components(full, ROTATIONAL_COMPONENTS, pair)
        singular = np.linalg.svd(subset, compute_uv=False)
        rank = int(np.sum(singular > max(subset.shape) * np.finfo(float).eps * (singular[0] if singular.size else 0.0)))
        cond = float(singular[0] / singular[-1]) if singular.size and singular[-1] > 0.0 else float("inf")
        score = (rank, cond)
        if score[0] > best_score[0] or (score[0] == best_score[0] and score[1] < best_score[1]):
            best = pair
            best_score = score
    return best


def _select_raw_components(raw: np.ndarray, component_names: tuple[str, ...], selected: tuple[str, ...]) -> np.ndarray:
    arr = np.asarray(raw, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape((-1, len(component_names)))
        idx = [component_names.index(item) for item in selected]
        return arr[:, idx].reshape(-1)
    idx = []
    for block in range(arr.shape[0] // len(component_names)):
        idx.extend(block * len(component_names) + component_names.index(item) for item in selected)
    return arr[idx, :]


def _constants_to_moments(constants: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(ROTCONST_TO_MOMENT / value if value > 0.0 else 0.0 for value in constants)


def _moment_weights(obs: IsotopologueObservation) -> tuple[float, float, float]:
    if obs.weights is None:
        return (1.0, 1.0, 1.0)
    constants = obs.corrected.as_tuple()
    sigmas_b = tuple((1.0 / weight) ** 0.5 for weight in obs.weights.as_tuple())
    weights = []
    for b_value, sigma_b in zip(constants, sigmas_b):
        sigma_i = abs(ROTCONST_TO_MOMENT * sigma_b / (b_value * b_value))
        weights.append(1.0 / (sigma_i * sigma_i) if sigma_i > 0.0 else 1.0)
    return tuple(weights)


def _predicate_observations(
    predicates: tuple[QMParameterPredicate, ...],
    labels: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray, list[tuple[str, str]]]:
    values = []
    weights = []
    row_labels = []
    for predicate in predicates:
        matches = _predicate_indices(predicate, labels)
        if not matches:
            raise ScientificValidationError(f"QM predicate did not match any GIC: {predicate.label_pattern}")
        for idx in matches:
            values.append(predicate.value)
            weights.append(predicate.weight)
            row_labels.append((predicate.source, labels[idx]))
    return np.array(values, dtype=float), np.array(weights, dtype=float), row_labels


def _predicate_values(
    predicates: tuple[QMParameterPredicate, ...],
    labels: tuple[str, ...],
    q_values: np.ndarray,
) -> np.ndarray:
    values = []
    for predicate in predicates:
        for idx in _predicate_indices(predicate, labels):
            values.append(float(q_values[idx]))
    return np.array(values, dtype=float)


def _predicate_indices(predicate: QMParameterPredicate, labels: tuple[str, ...]) -> list[int]:
    pattern = predicate.label_pattern.lower()
    return [idx for idx, label in enumerate(labels) if pattern in label.lower()]


def _is_planar(coords: np.ndarray, tol: float = 1.0e-3) -> bool:
    centered = np.asarray(coords, dtype=float) - np.mean(coords, axis=0)
    if centered.shape[0] < 4:
        return False
    singular = np.linalg.svd(centered, compute_uv=False)
    scale = max(float(singular[0]), 1.0)
    return float(singular[-1]) / scale < tol


def _residual_rows(
    model: MeasurementModel,
    calculated: np.ndarray,
    observed: np.ndarray,
) -> tuple[SemiexperimentalResidual, ...]:
    rows = []
    for idx, (isotopologue, label) in enumerate(model.labels):
        rows.append(SemiexperimentalResidual(isotopologue, label, float(observed[idx]), float(calculated[idx]), float(observed[idx] - calculated[idx])))
    return tuple(rows)


def _rotational_constant_rows(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
) -> tuple[SemiexperimentalRotationalConstantComparison, ...]:
    calculated = _constants_vector(atoms, coords, observations).reshape((-1, len(ROTATIONAL_COMPONENTS)))
    rows: list[SemiexperimentalRotationalConstantComparison] = []
    for obs, calc_triplet in zip(observations, calculated):
        for component, observed_value, calculated_value in zip(
            ROTATIONAL_COMPONENTS,
            obs.corrected.as_tuple(),
            calc_triplet,
        ):
            rows.append(
                SemiexperimentalRotationalConstantComparison(
                    obs.label,
                    component,
                    float(observed_value),
                    float(calculated_value),
                    float(observed_value - calculated_value),
                )
            )
    return tuple(rows)


def _parameters(
    labels: tuple[str, ...],
    values: np.ndarray,
    active_mask: np.ndarray,
    transform: np.ndarray | None = None,
    covariance: np.ndarray | None = None,
    class_by_gic: tuple[str, ...] = (),
) -> tuple[SemiexperimentalParameter, ...]:
    params = []
    active_positions = {idx: pos for pos, idx in enumerate(np.where(active_mask)[0])}
    for idx, label in enumerate(labels):
        active = bool(active_mask[idx])
        parameter_class = class_by_gic[idx] if idx < len(class_by_gic) else ""
        sigma = 0.0
        if active:
            pos = active_positions[idx]
            if transform is not None and covariance is not None and transform.size and covariance.size:
                row = np.asarray(transform[pos, :], dtype=float)
                if covariance.shape == (row.size, row.size):
                    variance = float(row @ covariance @ row)
                    sigma = float(np.sqrt(max(variance, 0.0)))
            elif covariance is not None and covariance.size and pos < covariance.shape[0]:
                sigma = float(np.sqrt(max(float(covariance[pos, pos]), 0.0)))
        params.append(SemiexperimentalParameter(label, float(values[idx]), sigma, active, parameter_class))
    return tuple(params)


def _covariance(jac: np.ndarray, residual: np.ndarray) -> np.ndarray:
    if jac.size == 0:
        return np.zeros((0, 0), dtype=float)
    dof = max(jac.shape[0] - jac.shape[1], 1)
    sigma2 = float(residual @ residual) / dof
    return sigma2 * np.linalg.pinv(jac.T @ jac, rcond=1.0e-10)

def _diagnostics(
    weighted_jac: np.ndarray,
    weighted_residual: np.ndarray,
    *,
    convergence_reason: str,
    damping: float,
    accepted_steps: int,
    rejected_steps: int,
    max_iterations: int,
    n_optimized_parameters: int,
    observable: str,
    components: tuple[str, ...],
    planar: bool,
    auto_pruned_parameters: tuple[str, ...] = (),
    prune_condition_target: float = 0.0,
    gicforge_calls: int = 0,
    coordinate_model_reuse_steps: int = 0,
    trust_radius: float = 0.0,
    last_trust_ratio: float = 0.0,
    last_line_search_scale: float = 0.0,
    b_projector_analytic_refreshes: int = 0,
    b_projector_secant_updates: int = 0,
    b_projector_secant_rejections: int = 0,
    last_b_projector_secant_error: float = 0.0,
    parameter_scale_min: float = 1.0,
    parameter_scale_max: float = 1.0,
    coordinate_model: str = "gic",
) -> SemiexperimentalFitDiagnostics:
    conditioning = rank_condition(weighted_jac)
    incremental_rank = _incremental_column_rank(weighted_jac)
    obj = objective(weighted_residual)
    dof = max(weighted_residual.size - weighted_jac.shape[1], 1) if weighted_jac.ndim == 2 else 1
    return SemiexperimentalFitDiagnostics(
        convergence_reason=convergence_reason,
        objective=obj,
        weighted_rms=float(np.sqrt(np.mean(weighted_residual * weighted_residual))) if weighted_residual.size else 0.0,
        reduced_chi_square=float((weighted_residual @ weighted_residual) / dof) if weighted_residual.size else 0.0,
        rank=conditioning.rank,
        incremental_rank=incremental_rank,
        condition_number=conditioning.condition_number,
        damping=float(damping),
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
        max_iterations=int(max_iterations),
        n_optimized_parameters=int(n_optimized_parameters),
        observable=observable,
        components=components,
        planar=planar,
        auto_pruned_parameters=auto_pruned_parameters,
        prune_condition_target=float(prune_condition_target),
        gicforge_calls=int(gicforge_calls),
        coordinate_model_reuse_steps=int(coordinate_model_reuse_steps),
        trust_radius=float(trust_radius),
        last_trust_ratio=float(last_trust_ratio),
        last_line_search_scale=float(last_line_search_scale),
        b_projector_analytic_refreshes=int(b_projector_analytic_refreshes),
        b_projector_secant_updates=int(b_projector_secant_updates),
        b_projector_secant_rejections=int(b_projector_secant_rejections),
        last_b_projector_secant_error=float(last_b_projector_secant_error),
        parameter_scale_min=float(parameter_scale_min),
        parameter_scale_max=float(parameter_scale_max),
        coordinate_model=coordinate_model,
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
