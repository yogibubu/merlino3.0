from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class WorkflowSpec:
    workflow_id: str
    title: str
    category: str
    service: str
    description: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    backends: tuple[str, ...] = ("python",)
    default_backend: str = "python"
    status: str = "planned"


def default_workflows() -> list[WorkflowSpec]:
    """Return the Merlino4 top-level workflow map."""
    return [
        WorkflowSpec(
            workflow_id="molecule",
            title="Molecule",
            category="Structure",
            service="merlino_geometry",
            description="Structure input, isotopes, topology, rings and symmetry.",
            inputs=("XYZ/xyzin/SMILES from Python side", "isotope selections"),
            outputs=("normalized molecular model", "topology/ring/symmetry report"),
            backends=("python",),
            status="legacy GUI available",
        ),
        WorkflowSpec(
            workflow_id="gic",
            title="GIC Definition / B Matrix",
            category="Coordinates",
            service="merlino_gic + merlino_fortran",
            description="Automatic GIC definition from Cartesian geometry plus reusable B-matrix evaluation from a frozen coordinate schema.",
            inputs=("reference Cartesian geometry", "current Cartesian geometry for B evaluation", "frozen GIC definition JSON"),
            outputs=("gic_definition.json", "Gaussian-readable GIC block", "GIC report", "B matrix CSV", "GIC values CSV"),
            backends=("python", "fortran77"),
            default_backend="fortran77",
            status="two-utility service available",
        ),
        WorkflowSpec(
            workflow_id="gic_gf",
            title="GIC Definition / B Matrix / GF-PED",
            category="Vibrations",
            service="merlino_gic + merlino_gf",
            description="Frozen GIC definition, reusable B-matrix evaluation, Cartesian Hessian reading, internal Hessian transformation, optional Pulay scaling, Wilson GF, frequencies, normal modes and PED.",
            inputs=("frozen GIC definition JSON", "Cartesian Hessian FCHK adapter", "optional current Cartesian geometry", "optional Pulay diagonal scaling factors"),
            outputs=("GF/PED report", "frequencies CSV", "G matrix CSV", "internal Hessian CSV", "normal modes CSV", "PED CSV", "manifest"),
            backends=("python",),
            status="frozen-GIC GF service available",
        ),
        WorkflowSpec(
            workflow_id="dvr",
            title="DVR",
            category="Dynamics",
            service="merlino_dvr + merlino_fortran",
            description="Gaussian scan/path output to DVR levels and wavefunctions.",
            inputs=("Gaussian log or grid CSV", "solver/boundary settings"),
            outputs=("levels CSV", "vectors/profile CSV", "summary", "manifest"),
            backends=("python", "fortran77"),
            default_backend="fortran77",
            status="legacy window available",
        ),
        WorkflowSpec(
            workflow_id="vpt2_vci",
            title="VPT2 / VCI",
            category="Vibrations",
            service="merlino_vpt2_vci",
            description="Anharmonic VPT2/VCI on canonical Merlino quartic force fields in Cartesian normal modes.",
            inputs=("canonical normal-mode QFF", "optional Gaussian FCHK frequency/QFF adapter", "basis cutoffs", "root/convergence settings"),
            outputs=("VPT2 levels", "VCI levels", "VPT2/VCI comparison", "dominant coefficients", "Davidson report"),
            backends=("python", "fortran77"),
            status="normal-mode anharmonic window available",
        ),
        WorkflowSpec(
            workflow_id="semiexp_geometry",
            title="Semiexperimental Geometry",
            category="Structure",
            service="merlino_semiexp",
            description="Cartesian equilibrium-geometry fit using automatically generated totally symmetric GICs or Hessian-free symmetry-adapted Cartesian displacements.",
            inputs=("parent Cartesian geometry", "isotopologue B0 constants", "QM vibrational/electronic corrections", "experimental uncertainties", "QM predicates", "shared/fixed parameter classes"),
            outputs=("fitted structure", "working-coordinate parameters", "primitive internal coordinates/errors", "residuals", "Kraitchman comparison", "covariance/correlation", "least-squares Hessian", "minimum check", "SVD/constraint diagnostics", "checkpoint/restart", "optional leave-one-out", "manifest"),
            backends=("python", "fortran77"),
            status="standard solver",
        ),
        WorkflowSpec(
            workflow_id="jobs_reports",
            title="Jobs / Reports",
            category="Project",
            service="merlino_core",
            description="Manifest browser, logs, reproducibility metadata and output collection.",
            inputs=("workflow manifests", "backend logs"),
            outputs=("status dashboard", "report bundle"),
            backends=("python",),
            status="planned",
        ),
    ]
