from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .gaussian_qff import anharmonic_input_from_gaussian_fchk, hessian_input_from_gaussian_fchk, read_indexed_qff_text
from .internal_gf import InternalGFResult, gf_from_hessian_input_with_merlino_gics
from .vci import QuarticForceField, VCIOptions, force_field_from_anharmonic_input
from .vpt2 import VPT2VCIComparison, compare_vpt2_vci


@dataclass(frozen=True)
class GFReport:
    fchk_path: Path
    result: InternalGFResult
    text: str


@dataclass(frozen=True)
class VPT2VCIReport:
    force_field: QuarticForceField
    comparison: VPT2VCIComparison
    text: str


def run_gf_report_from_fchk(fchk_path: Path) -> GFReport:
    """Read an FCHK adapter and return a formatted GF/PED report."""
    path = Path(fchk_path)
    hessian_input = hessian_input_from_gaussian_fchk(path)
    result = gf_from_hessian_input_with_merlino_gics(hessian_input)
    return GFReport(path, result, format_gf_report(path, result))


def load_force_field(fchk_path: Path | None = None, qff_path: Path | None = None) -> QuarticForceField:
    """Load a canonical Merlino QFF using optional Gaussian/FCHK frequencies."""
    frequencies = None
    anharmonic_input = None
    if fchk_path is not None:
        anharmonic_input = anharmonic_input_from_gaussian_fchk(Path(fchk_path))
        frequencies = (
            anharmonic_input.anharmonic_frequencies_cm
            if anharmonic_input.anharmonic_frequencies_cm.size
            else anharmonic_input.harmonic_frequencies_cm
        )
    if qff_path is not None:
        return read_indexed_qff_text(Path(qff_path), frequencies)
    if anharmonic_input is not None:
        return force_field_from_anharmonic_input(anharmonic_input)
    raise FileNotFoundError("Provide an existing FCHK file or indexed QFF text file")


def run_vpt2_vci_report(
    force_field: QuarticForceField,
    *,
    max_quanta: int,
    roots: int,
    options: VCIOptions | None = None,
) -> VPT2VCIReport:
    """Run VPT2/VCI on a canonical force field and return a formatted report."""
    comparison = compare_vpt2_vci(force_field, max_quanta=max_quanta, n_roots=roots, options=options)
    return VPT2VCIReport(force_field, comparison, format_vpt2_vci_report(force_field, comparison))


def format_gf_report(fchk_path: Path, result: InternalGFResult) -> str:
    lines = [
        "GF/PED from Merlino non-redundant GICs",
        f"Source FCHK: {Path(fchk_path)}",
        f"GIC count: {len(result.gic_labels)}",
        "",
        "Frequencies (cm-1):",
    ]
    for idx, freq in enumerate(result.frequencies_cm, start=1):
        lines.append(f"  mode {idx:3d}: {freq:12.3f}")

    lines.extend(["", "GIC labels:"])
    for idx, label in enumerate(result.gic_labels, start=1):
        lines.append(f"  GIC{idx:03d}: {label}")

    lines.extend(["", "PED (%) rows=GIC cols=modes:"])
    header = "          " + " ".join(f"M{idx:02d}" for idx in range(1, len(result.frequencies_cm) + 1))
    lines.append(header)
    for idx, row in enumerate(result.ped.values, start=1):
        values = " ".join(f"{value:7.2f}" for value in row)
        lines.append(f"  GIC{idx:03d} {values}")
    return "\n".join(lines)


def format_vpt2_vci_report(qff: QuarticForceField, comparison: VPT2VCIComparison) -> str:
    lines = [
        "VPT2/VCI comparison on canonical Merlino QFF",
        f"Modes used in input force field: {len(qff.harmonic_frequencies_cm)}",
        f"Cubic terms: {len(qff.cubic_cm)}",
        f"Quartic terms: {len(qff.quartic_cm)}",
        f"VCI basis size: {len(comparison.vci.basis)}",
        "Input harmonic frequencies (cm-1): "
        + ", ".join(f"{value:.3f}" for value in qff.harmonic_frequencies_cm),
        "",
        "Root     VPT2 abs      VCI abs        d_abs     VPT2 exc      VCI exc        d_exc",
    ]
    n = min(
        len(comparison.vpt2.energies_cm),
        len(comparison.vci.energies_cm),
        len(comparison.energy_differences_cm),
    )
    for idx in range(n):
        lines.append(
            f"{idx + 1:4d} "
            f"{comparison.vpt2.energies_cm[idx]:12.4f} "
            f"{comparison.vci.energies_cm[idx]:12.4f} "
            f"{comparison.energy_differences_cm[idx]:10.4f} "
            f"{comparison.vpt2.excitation_energies_cm[idx]:12.4f} "
            f"{comparison.vci.excitation_energies_cm[idx]:12.4f} "
            f"{comparison.excitation_differences_cm[idx]:10.4f}"
        )

    if comparison.vci.blocks:
        lines.extend(["", "Symmetry blocks:"])
        for block in comparison.vci.blocks:
            lines.append(f"  {block.label}: states={len(block.basis_indices)} roots={block.n_roots}")

    if comparison.vci.state_contributions:
        lines.extend(["", "Dominant VCI contributions:"])
        for root, contribution in enumerate(comparison.vci.state_contributions[:n], start=1):
            pieces = [
                f"{state}:{coeff:+.3f}"
                for state, coeff in contribution.dominant_basis_states[:4]
            ]
            lines.append(
                f"  root {root:3d}: <n>={np.array2string(contribution.mode_quanta, precision=3)} "
                + ", ".join(pieces)
            )
    return "\n".join(lines)
