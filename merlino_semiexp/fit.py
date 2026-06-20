from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
from io import StringIO
import re
import tempfile

import numpy as np

from geometry.inertia import principal_moments
from geometry.physical_constants import Phy, get_physical_constants
from geometry.rotational import rotational_constants_MHz
from geometry.structure import Structure
from merlino_core import ScientificValidationError, build_run_manifest
from merlino_gic import run_gicforge
from merlino_core.numerics import damped_normal_step, limit_step, objective, rank_condition
from topology.elements import atomic_symbol
from merlino_fit.survibfit.modify_geom import read_xyz, write_xyz
from merlino_fit.survibfit.pipeline import b_matrix_analytic
from merlino_fit.survibfit.primitives import Primitive, eval_primitives

from .contracts import IsotopologueObservation, ParameterClassConstraint, QMParameterPredicate, SemiexperimentalFitRequest
from .kraitchman import KraitchmanComparison, KraitchmanSeedResult, kraitchman_comparison, kraitchman_seed_geometry


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
    observable: str
    components: tuple[str, ...]
    planar: bool


@dataclass(frozen=True)
class MeasurementModel:
    observable: str
    components: tuple[str, ...]
    labels: tuple[tuple[str, str], ...]
    observed: np.ndarray
    weights: np.ndarray
    planar: bool


@dataclass
class GICForgeSEBackend:
    atoms: tuple[str, ...]
    root: Path
    counter: int = 0
    last_workdir: Path | None = None
    point_group: str | None = None

    def model(self, coords: np.ndarray):
        self.counter += 1
        workdir = self.root / f"iter_{self.counter:04d}"
        workdir.mkdir(parents=True, exist_ok=True)
        _write_gicforge_se_inputs(workdir, self.atoms, coords)
        result = run_gicforge(workdir)
        point_group = _gicforge_point_group(workdir / "provout")
        if self.point_group is None:
            self.point_group = point_group
        elif point_group != self.point_group:
            raise ScientificValidationError(
                f"GICForge point group changed from {self.point_group} to {point_group} in {workdir}"
            )
        gauin = result.files.get("gauin.symm") or result.files.get("gauin")
        if gauin is None:
            raise ScientificValidationError(f"GICForge did not produce gauin in {workdir}")
        self.last_workdir = workdir
        return _gicforge_gic_model(gauin)


@dataclass(frozen=True)
class SemiexperimentalFitResult:
    atoms: tuple[str, ...]
    initial_coordinates_angstrom: np.ndarray
    final_coordinates_angstrom: np.ndarray
    parameters: tuple[SemiexperimentalParameter, ...]
    residuals: tuple[SemiexperimentalResidual, ...]
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
    max_iter: int = 20,
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
    gicforge_backend = _make_gicforge_backend(tuple(atoms), outdir)

    prims, u_matrix, labels = _gic_model(coords, z_numbers, request, gicforge_backend)
    measurement_model = _build_measurement_model(request, atoms, coords, prims, u_matrix, labels)
    active_mask = _active_mask(labels, request.fixed_parameters) & _gicforge_a1_mask(labels)
    loop_max_iter = max_iter if np.any(active_mask) else 0

    current_damping = float(damping)
    accepted_steps = 0
    rejected_steps = 0
    convergence_reason = "max_iter" if loop_max_iter else "no_active_totally_symmetric_parameters"
    previous_objective = None
    iteration = 0
    for iteration in range(1, loop_max_iter + 1):
        prims, u_matrix, labels = _gic_model(coords, z_numbers, request, gicforge_backend)
        active_mask = _active_mask(labels, request.fixed_parameters, request.parameter_classes)
        active_mask &= _gicforge_a1_mask(labels)
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
        )
        transform, _reduced_names, _class_by_gic = _parameter_class_transform(labels, active_mask, request.parameter_classes)
        jac = jac_gic @ transform
        if np.sqrt(np.mean(residual * residual)) < tolerance_MHz:
            convergence_reason = "rms_tolerance"
            break
        jac_weighted = jac * sqrt_weights[:, None]
        gradient = jac_weighted.T @ weighted_residual
        if float(np.linalg.norm(gradient, ord=np.inf)) < gradient_tolerance:
            convergence_reason = "gradient_tolerance"
            break
        dq_reduced = damped_normal_step(jac_weighted, weighted_residual, current_damping)
        dq_reduced = limit_step(dq_reduced, max_step)
        dq_active = transform @ dq_reduced
        dq = np.zeros_like(q)
        dq[np.where(active_mask)[0]] = dq_active
        candidate, candidate_objective = _line_search_update(
            atoms,
            coords,
            z_numbers,
            request,
            labels,
            measurement_model,
            prims,
            u_matrix,
            dq,
            gicforge_backend,
            current_objective=current_objective,
        )
        if candidate_objective < current_objective:
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
        iteration = loop_max_iter

    prims, u_matrix, labels = _gic_model(coords, z_numbers, request, gicforge_backend)
    active_mask = _active_mask(labels, request.fixed_parameters, request.parameter_classes)
    active_mask &= _gicforge_a1_mask(labels)
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
        observable=measurement_model.observable,
        components=measurement_model.components,
        planar=measurement_model.planar,
    )
    sigmas_active = np.sqrt(np.clip(np.diag(covariance), 0.0, None)) if covariance.size else np.array(())
    parameters = _parameters(labels, q_final, active_mask, sigmas_active, transform, class_by_gic)
    residual_rows = _residual_rows(measurement_model, calc, obs)
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
            kraitchman_seed=kraitchman_seed,
            effective_parameter_names=reduced_names,
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


def write_semiexperimental_outputs(
    outdir: Path,
    request: SemiexperimentalFitRequest,
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    parameters: tuple[SemiexperimentalParameter, ...],
    residuals: tuple[SemiexperimentalResidual, ...],
    kraitchman: tuple[KraitchmanComparison, ...] = (),
    kraitchman_seed: KraitchmanSeedResult | None = None,
    effective_parameter_names: tuple[str, ...] = (),
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
    kraitchman_csv = outdir / "semiexp_kraitchman.csv"
    kraitchman_xyz = outdir / "semiexp_kraitchman_geometry.xyz"
    covariance_csv = outdir / "semiexp_covariance.csv"
    correlation_csv = outdir / "semiexp_correlation.csv"
    hessian_csv = outdir / "semiexp_hessian.csv"
    hessian_eigs_csv = outdir / "semiexp_hessian_eigenvalues.csv"
    diagnostics_csv = outdir / "semiexp_diagnostics.csv"
    active_names = effective_parameter_names or _effective_parameter_names(parameters)
    write_xyz(xyz, atoms, coords, comment="Merlino semiexperimental equilibrium geometry")
    params.write_text(parameters_csv(parameters), encoding="utf-8")
    residual_csv.write_text(residuals_csv(residuals), encoding="utf-8")
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
    manifest_inputs = {"initial_geometry": request.initial_geometry}
    coordinate_generation = {
        "primitive_source": "GICForge ReadAllGIC generated at each SE iteration",
        "reduction": "non-redundant GIC transform",
        "symmetry": "GICForge/symm.f same-type coordinate symmetrization with strict/quasi tolerance",
        "active_subspace": "GICForge-assigned A1 coordinates only",
        "ring_coordinates": "GICForge ring deformation and puckering coordinates",
        "gicforge_iterations": str(outdir / "gicforge_iterations"),
    }
    outputs = {
        "geometry": xyz,
        "parameters": params,
        "residuals": residual_csv,
        "kraitchman": kraitchman_csv,
        "covariance": covariance_csv,
        "correlation": correlation_csv,
        "hessian": hessian_csv,
        "hessian_eigenvalues": hessian_eigs_csv,
        "diagnostics": diagnostics_csv,
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
            "fixed_parameters": request.fixed_parameters,
            "parameter_classes": tuple(
                {"name": item.name, "patterns": item.patterns, "mode": item.mode}
                for item in request.parameter_classes
            ),
            "stationary_point": stationary_point,
            "convergence_reason": diagnostics.convergence_reason if diagnostics else "not_reported",
            "observable": diagnostics.observable if diagnostics else request.observable,
            "rotational_components": diagnostics.components if diagnostics else request.rotational_components,
            "isotopologues": tuple(obs.label for obs in request.observations),
            "n_isotopologues": len(request.observations),
            "n_qm_predicates": len(request.qm_predicates),
            "n_gic_parameters": len(parameters),
            "n_effective_parameters": len(active_names),
            "n_active_gic_parameters": sum(1 for item in parameters if item.active),
            "n_kraitchman_rows": len(kraitchman),
            "kraitchman_seed_method": kraitchman_seed.method if kraitchman_seed else "not_available",
            "n_kraitchman_seed_atoms": len(kraitchman_seed.fitted_atom_indices) if kraitchman_seed else 0,
            "rank": diagnostics.rank if diagnostics else None,
            "condition_number": diagnostics.condition_number if diagnostics else None,
            "weighted_rms": diagnostics.weighted_rms if diagnostics else None,
            "reduced_chi_square": diagnostics.reduced_chi_square if diagnostics else None,
            "coordinate_generation": coordinate_generation,
        },
        backend={
            "solver": "python-orchestrated",
            "coordinate_model": "gicforge-iterative-readallgic",
            "b_matrix": "analytic",
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


def _make_gicforge_backend(atoms: tuple[str, ...], outdir: Path | None) -> GICForgeSEBackend:
    if outdir is None:
        root = Path(tempfile.mkdtemp(prefix="merlino_se_gicforge_"))
    else:
        root = Path(outdir) / "gicforge_iterations"
        root.mkdir(parents=True, exist_ok=True)
    return GICForgeSEBackend(atoms=atoms, root=root)


def _write_gicforge_se_inputs(workdir: Path, atoms: tuple[str, ...], coords: np.ndarray) -> None:
    (workdir / "provin").write_text(
        "# GNIC SYMMALL BMAT ECKART G16 CLEAN\n\n"
        "Merlino semiexperimental GICForge coordinates\n\n"
        "0 1\n",
        encoding="utf-8",
    )
    lines = [str(len(atoms)), ""]
    for atom, (x, y, z) in zip(atoms, coords):
        lines.append(f"{atom:>4s} {x: 16.8E} {y: 16.8E} {z: 16.8E}")
    (workdir / "xyzin").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _gicforge_point_group(provout: Path) -> str:
    text = provout.read_text(encoding="utf-8", errors="replace") if provout.exists() else ""
    match = re.search(r"Point Group from symm\.f:\s*([A-Za-z0-9]+)", text)
    return match.group(1) if match else "UNKNOWN"


def _gicforge_gic_model(gauin: Path):
    irreps = _read_gicforge_irreps(Path(gauin).with_name("gicsym"))
    prims: list[Primitive] = []
    prim_index: dict[Primitive, int] = {}
    columns: list[np.ndarray] = []
    labels: list[str] = []
    for raw in Path(gauin).read_text(encoding="utf-8", errors="replace").splitlines():
        parsed = _parse_gicforge_line(raw)
        if parsed is None:
            continue
        name, terms, expression = parsed
        column = np.zeros(len(prims), dtype=float)
        for coeff, primitive in terms:
            if primitive not in prim_index:
                prim_index[primitive] = len(prims)
                prims.append(primitive)
                column = np.pad(column, (0, 1))
                for idx, existing in enumerate(columns):
                    columns[idx] = np.pad(existing, (0, 1))
            column[prim_index[primitive]] += coeff
        label_index = len(labels) + 1
        irrep = irreps.get(name, "UNK")
        labels.append(f"GIC{label_index:03d} GICForge {name} irrep={irrep} {expression} {_gicforge_aliases(terms)}")
        columns.append(column)
    if not columns:
        raise ScientificValidationError(f"No linear GICForge coordinates found in {gauin}")
    u_matrix = np.column_stack(columns)
    return prims, u_matrix, tuple(labels)


def _parse_gicforge_line(line: str) -> tuple[str, list[tuple[float, Primitive]], str] | None:
    stripped = line.strip()
    if not stripped or "=" not in stripped:
        return None
    name = stripped.split("=", 1)[0].strip()
    rhs = stripped.split("=", 1)[1].strip()
    terms: list[tuple[float, Primitive]] = []
    number = r"[+-]?\s*(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?"
    for match in re.finditer(rf"({number})\s*\*\s*([RADLU])\(([^)]*)\)", rhs):
        coeff = float(match.group(1).replace(" ", "").replace("D", "E").replace("d", "e"))
        primitive = _gicforge_primitive(match.group(2), match.group(3))
        terms.append((coeff, primitive))
    if not terms:
        simple = re.search(r"\b([RADLU])\(([^)]*)\)", rhs)
        if simple:
            terms.append((1.0, _gicforge_primitive(simple.group(1), simple.group(2))))
    if not terms:
        return None
    return name.replace("(Inactive)", "").strip(), terms, rhs


def _read_gicforge_irreps(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    irreps: dict[str, str] = {}
    for idx, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines()):
        if idx == 0 and raw.lower().startswith("name,"):
            continue
        parts = [part.strip() for part in raw.split(",")]
        if len(parts) >= 2 and parts[0]:
            irreps[parts[0]] = parts[1]
    return irreps


def _gicforge_primitive(kind: str, atoms_text: str) -> Primitive:
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
        mode_token = values[4]
        mode = mode_token if mode_token in {-1, -2} else -1
        return Primitive("linear_bend", atoms[:3], mode=mode)
    raise ScientificValidationError(f"Unsupported GICForge primitive {kind}({atoms_text})")


def _gicforge_aliases(terms: list[tuple[float, Primitive]]) -> str:
    aliases = []
    for _coeff, primitive in terms:
        aliases.append(_primitive_alias(primitive))
    return " ".join(sorted(set(aliases)))


def _primitive_alias(primitive: Primitive) -> str:
    atoms = ",".join(str(atom + 1) for atom in primitive.atoms)
    if primitive.kind == "bond":
        return f"bond({atoms})"
    if primitive.kind == "angle":
        return f"angle({atoms})"
    if primitive.kind == "dihedral":
        return f"dihedral({atoms})"
    if primitive.kind == "out_of_plane":
        return f"out_of_plane({atoms})"
    if primitive.kind == "linear_bend":
        return f"linear_bend({atoms};mode={primitive.mode})"
    return f"{primitive.kind}({atoms})"


def _gicforge_a1_mask(labels: tuple[str, ...]) -> np.ndarray:
    irreps = []
    for label in labels:
        match = re.search(r"\birrep=([A-Za-z0-9'\"+-]+)", label)
        irreps.append(match.group(1) if match else None)
    if not any(irrep is not None for irrep in irreps):
        return np.ones(len(labels), dtype=bool)
    return np.array([irrep in {"A1", "A", "Ag", "A'"} for irrep in irreps], dtype=bool)


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
) -> np.ndarray:
    active_indices = np.where(active_mask)[0]
    base_q = _gic_values(prims, u_matrix, coords)
    jac = np.zeros((len(measurement_model.observed), len(active_indices)), dtype=float)
    for col, idx in enumerate(active_indices):
        dq = np.zeros_like(base_q)
        dq[idx] = step
        plus = _displace_along_gics(coords, prims, u_matrix, dq)
        plus_q = _gic_values(prims, u_matrix, plus)
        dq[idx] = -step
        minus = _displace_along_gics(coords, prims, u_matrix, dq)
        minus_q = _gic_values(prims, u_matrix, minus)
        jac[:, col] = (
            _measurement_vector(atoms, plus, request, plus_q, labels, measurement_model)
            - _measurement_vector(atoms, minus, request, minus_q, labels, measurement_model)
        ) / (2.0 * step)
    return jac


def _displace_along_gics(coords: np.ndarray, prims: object, u_matrix: np.ndarray, dq: np.ndarray) -> np.ndarray:
    bq = u_matrix.T @ b_matrix_analytic(prims, coords)
    dx = np.linalg.pinv(bq, rcond=1.0e-8) @ dq
    return coords + dx.reshape(coords.shape)


def _line_search_update(
    atoms: list[str],
    coords: np.ndarray,
    z_numbers: np.ndarray,
    request: SemiexperimentalFitRequest,
    labels: tuple[str, ...],
    measurement_model: "MeasurementModel",
    prims: object,
    u_matrix: np.ndarray,
    dq: np.ndarray,
    gicforge_backend: GICForgeSEBackend,
    *,
    current_objective: float,
) -> tuple[np.ndarray, float]:
    observed = measurement_model.observed
    sqrt_weights = np.sqrt(measurement_model.weights)
    best_coords = coords
    best_objective = current_objective
    for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625):
        candidate = _displace_along_gics(coords, prims, u_matrix, scale * dq)
        try:
            cand_prims, cand_u_matrix, cand_labels = _gic_model(candidate, z_numbers, request, gicforge_backend)
            q_candidate = _gic_values(cand_prims, cand_u_matrix, candidate)
            calc = _measurement_vector(atoms, candidate, request, q_candidate, cand_labels, measurement_model)
        except Exception:
            continue
        if calc.shape != observed.shape:
            continue
        residual = (observed - calc) * sqrt_weights
        candidate_objective = objective(residual)
        if candidate_objective < best_objective:
            best_coords = candidate
            best_objective = candidate_objective
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
    base_q = _gic_values(prims, u_matrix, coords)
    full = np.zeros((3 * len(observations), len(base_q)), dtype=float)
    for idx in range(len(base_q)):
        dq = np.zeros_like(base_q)
        dq[idx] = 1.0e-4
        plus = _displace_along_gics(coords, prims, u_matrix, dq)
        dq[idx] = -1.0e-4
        minus = _displace_along_gics(coords, prims, u_matrix, dq)
        full[:, idx] = (_constants_vector(atoms, plus, observations) - _constants_vector(atoms, minus, observations)) / 2.0e-4
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


def _parameters(
    labels: tuple[str, ...],
    values: np.ndarray,
    active_mask: np.ndarray,
    sigmas_active: np.ndarray,
    transform: np.ndarray | None = None,
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
            if transform is not None and transform.size:
                cols = np.where(np.abs(transform[pos, :]) > 0.0)[0]
                if cols.size and cols[0] < len(sigmas_active):
                    sigma = float(abs(transform[pos, cols[0]]) * sigmas_active[cols[0]])
            elif pos < len(sigmas_active):
                sigma = float(sigmas_active[pos])
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
    observable: str,
    components: tuple[str, ...],
    planar: bool,
) -> SemiexperimentalFitDiagnostics:
    conditioning = rank_condition(weighted_jac)
    obj = objective(weighted_residual)
    dof = max(weighted_residual.size - weighted_jac.shape[1], 1) if weighted_jac.ndim == 2 else 1
    return SemiexperimentalFitDiagnostics(
        convergence_reason=convergence_reason,
        objective=obj,
        weighted_rms=float(np.sqrt(np.mean(weighted_residual * weighted_residual))) if weighted_residual.size else 0.0,
        reduced_chi_square=float((weighted_residual @ weighted_residual) / dof) if weighted_residual.size else 0.0,
        rank=conditioning.rank,
        condition_number=conditioning.condition_number,
        damping=float(damping),
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
        observable=observable,
        components=components,
        planar=planar,
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
