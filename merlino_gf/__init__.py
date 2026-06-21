"""Wilson GF/PED services for Merlino."""

from importlib import import_module

from .harmonic import GFResult, mass_weighted_cartesian_hessian, solve_wilson_gf
from .internal import (
    BOHR_TO_ANGSTROM,
    InternalGFResult,
    PEDTable,
    gf_from_cartesian_hessian_and_gic_b_matrix,
    gf_from_cartesian_hessian_and_merlino_gics,
    gf_from_gaussian_fchk_with_merlino_gics,
    gf_from_hessian_input_and_gic_definition,
    gf_from_hessian_input_with_merlino_gics,
    gic_labels_from_u,
    primitive_label,
    pulay_scale_internal_hessian,
)
from .models import HessianInput
_SERVICE_EXPORTS = {
    "GFReport",
    "format_gf_report",
    "gf_csv_tables",
    "pulay_scaling_factors",
    "run_gf_report_from_fchk",
    "run_gic_gf_report_from_fchk",
    "write_csv_tables",
}


def __getattr__(name: str):
    if name in _SERVICE_EXPORTS:
        return getattr(import_module(".service", __name__), name)
    raise AttributeError(name)

__all__ = [
    "BOHR_TO_ANGSTROM",
    "GFReport",
    "GFResult",
    "HessianInput",
    "InternalGFResult",
    "PEDTable",
    "format_gf_report",
    "gf_csv_tables",
    "gf_from_cartesian_hessian_and_gic_b_matrix",
    "gf_from_cartesian_hessian_and_merlino_gics",
    "gf_from_gaussian_fchk_with_merlino_gics",
    "gf_from_hessian_input_and_gic_definition",
    "gf_from_hessian_input_with_merlino_gics",
    "gic_labels_from_u",
    "mass_weighted_cartesian_hessian",
    "primitive_label",
    "pulay_scale_internal_hessian",
    "pulay_scaling_factors",
    "run_gf_report_from_fchk",
    "run_gic_gf_report_from_fchk",
    "solve_wilson_gf",
    "write_csv_tables",
]
