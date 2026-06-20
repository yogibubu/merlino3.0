from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

import numpy as np

from merlino_fit.survibfit.modify_geom import read_xyz

from .contracts import IsotopologueObservation, ParameterClassConstraint, SemiexperimentalFitRequest
from .fit import SemiexperimentalFitResult, _atomic_number, _gic_model, fit_semiexperimental_geometry


@dataclass(frozen=True)
class SemiexperimentalGICPreview:
    atoms: tuple[str, ...]
    gic_labels: tuple[str, ...]
    suggested_classes: tuple[ParameterClassConstraint, ...]
    warnings: tuple[str, ...]

    @property
    def text(self) -> str:
        lines = [
            "Merlino semiexperimental GIC preview",
            f"Atoms: {len(self.atoms)}",
            f"Non-redundant GICs: {len(self.gic_labels)}",
            "",
            "Suggested parameter classes:",
        ]
        lines.extend(
            f"  {item.name}:{item.mode}:{'|'.join(item.patterns)}" for item in self.suggested_classes
        )
        if not self.suggested_classes:
            lines.append("  none")
        lines.extend(["", "GIC labels:"])
        lines.extend(f"  {label}" for label in self.gic_labels)
        if self.warnings:
            lines.extend(["", "Warnings:", *[f"  {item}" for item in self.warnings]])
        return "\n".join(lines)


@dataclass(frozen=True)
class SemiexperimentalBenchmarkCase:
    label: str
    request: SemiexperimentalFitRequest


@dataclass(frozen=True)
class SemiexperimentalBenchmarkRow:
    label: str
    rms_MHz: float
    iterations: int
    rank: int
    condition_number: float
    stationary_point: str
    n_parameters: int
    n_kraitchman: int


def preview_semiexperimental_gics(
    xyz: Path,
    observations: tuple[IsotopologueObservation, ...] = (),
) -> SemiexperimentalGICPreview:
    atoms, coords, _comment = read_xyz(Path(xyz))
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    _prims, _u_matrix, labels = _gic_model(np.asarray(coords, dtype=float), z_numbers)
    suggestions = suggest_parameter_classes(tuple(atoms), labels, observations)
    warnings = _preview_warnings(labels, suggestions)
    return SemiexperimentalGICPreview(tuple(atoms), labels, suggestions, warnings)


def suggest_parameter_classes(
    atoms: tuple[str, ...],
    gic_labels: tuple[str, ...],
    observations: tuple[IsotopologueObservation, ...] = (),
) -> tuple[ParameterClassConstraint, ...]:
    suggestions: list[ParameterClassConstraint] = []
    h_substituted = _substituted_hydrogens(atoms, observations)
    h_atoms = tuple(idx + 1 for idx, atom in enumerate(atoms) if atom.upper() == "H")
    if h_atoms:
        for heavy in sorted({heavy for heavy, h in _heavy_h_bonds(atoms, gic_labels)}):
            patterns = tuple(f"bond({heavy},{h})" for h in h_atoms if _matches_any(gic_labels, f"bond({heavy},{h})"))
            if len(patterns) >= 2:
                name = f"{atoms[heavy - 1].upper()}H_stretches"
                suggestions.append(ParameterClassConstraint(name, patterns, "shared"))
        angle_patterns = tuple(label.split(maxsplit=1)[1] for label in gic_labels if "angle(" in label and _angle_has_h(label, atoms))
        if angle_patterns and not h_substituted:
            suggestions.append(ParameterClassConstraint("XH_angles", angle_patterns, "fixed"))
    return tuple(suggestions)


def write_semiexperimental_html_report(
    path: Path,
    result: SemiexperimentalFitResult,
    request: SemiexperimentalFitRequest,
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    html = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>Merlino semiexperimental geometry report</title>",
        "<style>body{font-family:Helvetica,Arial,sans-serif;margin:32px;line-height:1.35}"
        "table{border-collapse:collapse;margin:16px 0;width:100%}"
        "th,td{border:1px solid #ccc;padding:5px 7px;text-align:left;font-size:13px}"
        "th{background:#f3f3f3} code{background:#f7f7f7;padding:1px 3px}</style>",
        "</head><body>",
        "<h1>Merlino Semiexperimental Geometry Report</h1>",
        "<h2>Diagnostics</h2>",
        "<table><tr><th>Quantity</th><th>Value</th></tr>",
        _row("RMS", f"{result.rms_MHz:.8g}"),
        _row("Iterations", str(result.iterations)),
        _row("Stationary point", result.stationary_point),
        _row("Convergence", result.diagnostics.convergence_reason),
        _row("Rank", str(result.diagnostics.rank)),
        _row("Condition number", f"{result.diagnostics.condition_number:.8g}"),
        _row("Observable", result.diagnostics.observable),
        _row("Components", ",".join(result.diagnostics.components)),
        "</table>",
        "<h2>Parameter Classes</h2>",
        _classes_table(request.parameter_classes),
        "<h2>GIC Parameters</h2>",
        _parameters_table(result),
        "<h2>Residuals</h2>",
        _residuals_table(result),
        "<h2>Kraitchman Comparison</h2>",
        _kraitchman_table(result),
        "</body></html>",
    ]
    target.write_text("\n".join(html) + "\n", encoding="utf-8")
    return target


def run_semiexperimental_benchmark(
    cases: tuple[SemiexperimentalBenchmarkCase, ...],
    *,
    outdir: Path | None = None,
    max_iter: int = 12,
) -> tuple[SemiexperimentalBenchmarkRow, ...]:
    rows = []
    for case in cases:
        case_out = Path(outdir) / case.label if outdir is not None else None
        result = fit_semiexperimental_geometry(case.request, max_iter=max_iter, outdir=case_out)
        rows.append(
            SemiexperimentalBenchmarkRow(
                case.label,
                result.rms_MHz,
                result.iterations,
                result.diagnostics.rank,
                result.diagnostics.condition_number,
                result.stationary_point,
                len(result.parameters),
                len(result.kraitchman),
            )
        )
    return tuple(rows)


def benchmark_csv(rows: tuple[SemiexperimentalBenchmarkRow, ...]) -> str:
    lines = ["label,rms,iterations,rank,condition_number,stationary_point,n_parameters,n_kraitchman"]
    for row in rows:
        lines.append(
            f"{row.label},{row.rms_MHz:.12g},{row.iterations},{row.rank},"
            f"{row.condition_number:.12g},{row.stationary_point},{row.n_parameters},{row.n_kraitchman}"
        )
    return "\n".join(lines) + "\n"


def _preview_warnings(
    labels: tuple[str, ...],
    suggestions: tuple[ParameterClassConstraint, ...],
) -> tuple[str, ...]:
    warnings = []
    if not labels:
        warnings.append("No non-redundant GIC labels were generated")
    if not suggestions:
        warnings.append("No automatic parameter-class suggestion was found")
    return tuple(warnings)


def _substituted_hydrogens(atoms: tuple[str, ...], observations: tuple[IsotopologueObservation, ...]) -> set[int]:
    result = set()
    for obs in observations:
        for atom_index in obs.substitutions:
            if 1 <= atom_index <= len(atoms) and atoms[atom_index - 1].upper() == "H":
                result.add(atom_index)
    return result


def _heavy_h_bonds(atoms: tuple[str, ...], labels: tuple[str, ...]) -> tuple[tuple[int, int], ...]:
    pairs = []
    for label in labels:
        for left, right in _label_atom_pairs(label, "bond"):
            if atoms[left - 1].upper() == "H" and atoms[right - 1].upper() != "H":
                pairs.append((right, left))
            elif atoms[right - 1].upper() == "H" and atoms[left - 1].upper() != "H":
                pairs.append((left, right))
    return tuple(pairs)


def _label_atom_pairs(label: str, kind: str) -> tuple[tuple[int, int], ...]:
    pairs = []
    marker = f"{kind}("
    start = 0
    while True:
        pos = label.find(marker, start)
        if pos < 0:
            break
        end = label.find(")", pos)
        if end < 0:
            break
        parts = [int(part.strip()) for part in label[pos + len(marker):end].split(",")]
        if len(parts) >= 2:
            pairs.append((parts[0], parts[-1]))
        start = end + 1
    return tuple(pairs)


def _matches_any(labels: tuple[str, ...], pattern: str) -> bool:
    low = pattern.lower()
    return any(low in label.lower() for label in labels)


def _angle_has_h(label: str, atoms: tuple[str, ...]) -> bool:
    marker = "angle("
    pos = label.find(marker)
    if pos < 0:
        return False
    end = label.find(")", pos)
    if end < 0:
        return False
    atom_ids = [int(part.strip()) for part in label[pos + len(marker):end].split(",")]
    return any(1 <= idx <= len(atoms) and atoms[idx - 1].upper() == "H" for idx in atom_ids)


def _row(key: str, value: str) -> str:
    return f"<tr><th>{escape(key)}</th><td>{escape(value)}</td></tr>"


def _classes_table(classes: tuple[ParameterClassConstraint, ...]) -> str:
    if not classes:
        return "<p>No parameter classes.</p>"
    rows = ["<table><tr><th>Name</th><th>Mode</th><th>Patterns</th></tr>"]
    for item in classes:
        rows.append(
            f"<tr><td>{escape(item.name)}</td><td>{escape(item.mode)}</td>"
            f"<td><code>{escape('|'.join(item.patterns))}</code></td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _parameters_table(result: SemiexperimentalFitResult) -> str:
    rows = ["<table><tr><th>Name</th><th>Value</th><th>Sigma</th><th>Active</th><th>Class</th></tr>"]
    for item in result.parameters:
        rows.append(
            f"<tr><td>{escape(item.name)}</td><td>{item.value:.10g}</td><td>{item.sigma:.10g}</td>"
            f"<td>{int(item.active)}</td><td>{escape(item.parameter_class)}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _residuals_table(result: SemiexperimentalFitResult) -> str:
    rows = ["<table><tr><th>Isotopologue</th><th>Observable</th><th>Observed</th><th>Calculated</th><th>Residual</th></tr>"]
    for item in result.residuals:
        rows.append(
            f"<tr><td>{escape(item.isotopologue)}</td><td>{escape(item.constant)}</td>"
            f"<td>{item.observed_equilibrium_MHz:.10g}</td><td>{item.calculated_MHz:.10g}</td>"
            f"<td>{item.residual_MHz:.10g}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _kraitchman_table(result: SemiexperimentalFitResult) -> str:
    if not result.kraitchman:
        return "<p>No single-substitution Kraitchman comparison available.</p>"
    rows = ["<table><tr><th>Isotopologue</th><th>Atom</th><th>Axis</th><th>Kraitchman abs A</th><th>Fit abs A</th><th>Difference A</th></tr>"]
    for item in result.kraitchman:
        rows.append(
            f"<tr><td>{escape(item.isotopologue)}</td><td>{item.atom_index} {escape(item.atom)}</td>"
            f"<td>{escape(item.coordinate)}</td><td>{item.kraitchman_abs_angstrom:.10g}</td>"
            f"<td>{item.fitted_abs_angstrom:.10g}</td><td>{item.difference_angstrom:.10g}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)
