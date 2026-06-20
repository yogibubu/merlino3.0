from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class WorkflowSpec:
    workflow_id: str
    title: str
    service: str
    description: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    status: str = "planned"


def default_workflows() -> list[WorkflowSpec]:
    """Return the Merlino4 top-level workflow map."""
    return [
        WorkflowSpec(
            workflow_id="molecule",
            title="Molecule",
            service="merlino_geometry",
            description="Structure input, isotopes, topology, rings and symmetry.",
            inputs=("XYZ/xyzin/SMILES from Python side", "isotope selections"),
            outputs=("normalized molecular model", "topology/ring/symmetry report"),
            status="legacy GUI available",
        ),
        WorkflowSpec(
            workflow_id="gic",
            title="GIC / Gaussian Input",
            service="merlino_gic + merlino_fortran",
            description="GIC construction, GICForge execution and Gaussian GIC input.",
            inputs=("cartesian XYZ", "GIC options", "Gaussian route/resources"),
            outputs=("gauin.gjf", "GIC report", "optional B matrix", "manifest"),
            status="service extraction in progress",
        ),
        WorkflowSpec(
            workflow_id="dvr",
            title="DVR",
            service="merlino_dvr + merlino_fortran",
            description="Gaussian scan/path output to DVR levels and wavefunctions.",
            inputs=("Gaussian log or grid CSV", "solver/boundary settings"),
            outputs=("levels CSV", "vectors/profile CSV", "summary", "manifest"),
            status="legacy window available",
        ),
        WorkflowSpec(
            workflow_id="vpt2_vci",
            title="VPT2 / VCI",
            service="merlino_vpt2_vci",
            description="VPT2/VCI from Gaussian quartic force fields with Davidson for large VCI spaces.",
            inputs=("Gaussian quartic force field", "basis cutoffs", "root/convergence settings"),
            outputs=("VPT2 constants", "VCI levels", "dominant coefficients", "Davidson report"),
        ),
        WorkflowSpec(
            workflow_id="semiexp_geometry",
            title="Semiexperimental Geometry",
            service="merlino_semiexp",
            description="Equilibrium geometry from isotopologue rotational constants and QM vibrational corrections.",
            inputs=("experimental rotational constants", "isotopologues", "QM vibrational corrections"),
            outputs=("fitted structure", "residuals", "covariance/correlation", "manifest"),
        ),
        WorkflowSpec(
            workflow_id="jobs_reports",
            title="Jobs / Reports",
            service="merlino_core",
            description="Manifest browser, logs, reproducibility metadata and output collection.",
            inputs=("workflow manifests", "backend logs"),
            outputs=("status dashboard", "report bundle"),
            status="planned",
        ),
    ]
