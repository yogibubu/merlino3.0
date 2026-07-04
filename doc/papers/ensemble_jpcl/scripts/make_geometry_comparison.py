from __future__ import annotations

import csv
import json
import math
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from merlino_semiexp.ensemble import (
    _atomic_number,
    _class_projection,
    _ensemble_classes_from_mapping,
    _ensemble_molecules_from_mapping,
    _normalize_kind,
    _primitive_terms,
    _synthon_signatures,
    _synthon_zeff,
)
from merlino_semiexp.geometry_input import read_geometry_input


PAPER = ROOT / "doc" / "papers" / "ensemble_jpcl"
ANALYSIS = PAPER / "analysis" / "geometry_comparison"
GENERATED = PAPER / "generated"


@dataclass(frozen=True)
class Row:
    set_name: str
    molecule: str
    class_name: str
    kind: str
    label: str
    qc: float
    ensemble: float
    separate: float
    unit: str

    @property
    def delta(self) -> float:
        return self.ensemble - self.separate


def _load_corrections(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["class"]: float(row["correction"]) for row in csv.DictReader(handle)}


def _primitive_value(label: str, coords: np.ndarray) -> tuple[str, float] | None:
    terms = _primitive_terms(label)
    if len(terms) != 1:
        return None
    kind, indices, coeff = terms[0]
    if abs(coeff - 1.0) > 1.0e-12:
        return None
    if kind == "stretch" and len(indices) == 2:
        i, j = indices[0] - 1, indices[1] - 1
        return kind, float(np.linalg.norm(coords[i] - coords[j]))
    if kind == "bend" and len(indices) == 3:
        i, j, k = indices[0] - 1, indices[1] - 1, indices[2] - 1
        v1 = coords[i] - coords[j]
        v2 = coords[k] - coords[j]
        denom = np.linalg.norm(v1) * np.linalg.norm(v2)
        if denom <= 0.0:
            return None
        cosine = float(np.clip(np.dot(v1, v2) / denom, -1.0, 1.0))
        return kind, math.degrees(math.acos(cosine))
    return None


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _geometry_rows(
    *,
    set_name: str,
    ensemble_job: Path,
    corrections_csv: Path,
    separate_dirs: dict[str, Path],
) -> list[Row]:
    data = tomllib.loads(ensemble_job.read_text(encoding="utf-8"))
    molecules = _ensemble_molecules_from_mapping(ensemble_job, data["molecules"])
    classes = _ensemble_classes_from_mapping(data["classes"])
    corrections = _load_corrections(corrections_csv)
    rows: list[Row] = []

    for molecule in molecules:
        if molecule.name not in separate_dirs:
            raise KeyError(f"Missing separate-fit directory for {molecule.name}")
        geometry = read_geometry_input(Path(molecule.request.initial_geometry))
        atoms = tuple(geometry.atoms)
        coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
        z_numbers = np.asarray([_atomic_number(atom) for atom in atoms], dtype=int)
        synthon_signatures = _synthon_signatures(coords, z_numbers)
        synthon_values = _synthon_zeff(coords, z_numbers)
        parameter_csv = separate_dirs[molecule.name] / "semiexp_geometry_parameters.csv"
        if not parameter_csv.exists():
            raise FileNotFoundError(parameter_csv)
        for param in _csv_rows(parameter_csv):
            primitive = _primitive_value(param["label"], coords)
            if primitive is None:
                continue
            primitive_kind, qc_value = primitive
            if primitive_kind == "stretch":
                separate_text = param.get("value_angstrom", "")
                unit = "A"
            elif primitive_kind == "bend":
                separate_text = param.get("value_degree", "")
                unit = "deg"
            else:
                continue
            if not separate_text:
                continue
            matching = []
            for item in classes:
                if _normalize_kind(item.kind) != primitive_kind:
                    continue
                projection = _class_projection(
                    item,
                    param["label"],
                    atoms,
                    coords=coords,
                    synthon_signatures=synthon_signatures,
                    synthon_zeff=synthon_values,
                )
                if abs(projection) > 1.0e-12:
                    matching.append((item.name, projection))
            if len(matching) != 1:
                continue
            class_name, projection = matching[0]
            correction = projection * corrections[class_name]
            if unit == "deg":
                correction = math.degrees(correction)
            rows.append(
                Row(
                    set_name=set_name,
                    molecule=molecule.name,
                    class_name=class_name,
                    kind=primitive_kind,
                    label=param["label"],
                    qc=qc_value,
                    ensemble=qc_value + correction,
                    separate=float(separate_text),
                    unit=unit,
                )
            )
    return rows


def _write_detail_csv(path: Path, rows: list[Row]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["set", "molecule", "class", "kind", "label", "qc", "ensemble", "separate", "delta", "unit"])
        for row in rows:
            writer.writerow(
                [
                    row.set_name,
                    row.molecule,
                    row.class_name,
                    row.kind,
                    row.label,
                    f"{row.qc:.10f}",
                    f"{row.ensemble:.10f}",
                    f"{row.separate:.10f}",
                    f"{row.delta:.10f}",
                    row.unit,
                ]
            )


def _aggregate(rows: list[Row]) -> list[tuple[str, str, str, int, float, float, float, float, float]]:
    grouped: dict[tuple[str, str, str], list[Row]] = {}
    for row in rows:
        grouped.setdefault((row.class_name, row.kind, row.unit), []).append(row)
    output = []
    for (class_name, kind, unit), items in sorted(grouped.items(), key=lambda item: (item[0][1], item[0][0])):
        qcs = np.asarray([item.qc for item in items], dtype=float)
        ensembles = np.asarray([item.ensemble for item in items], dtype=float)
        separates = np.asarray([item.separate for item in items], dtype=float)
        deltas = ensembles - separates
        output.append(
            (
                class_name.replace("_", r"\_"),
                kind,
                unit,
                len(items),
                float(np.mean(qcs)),
                float(np.mean(ensembles)),
                float(np.mean(separates)),
                float(np.mean(deltas)),
                float(np.max(np.abs(deltas))),
            )
        )
    return output


def _write_table(path: Path, rows: list[Row], *, title: str) -> None:
    aggregate = _aggregate(rows)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(r"\begin{tabular}{llrrrrr}" + "\n")
        handle.write(r"\toprule" + "\n")
        handle.write(r"Class & Unit & $n$ & QC & ensemble & separate & $\Delta_\mathrm{ens-sep}$ \\" + "\n")
        handle.write(r"\midrule" + "\n")
        for class_name, _kind, unit, n_items, qc, ensemble, separate, delta, max_abs in aggregate:
            if unit == "A":
                qc_text = f"{qc:.4f}"
                ensemble_text = f"{ensemble:.4f}"
                separate_text = f"{separate:.4f}"
                delta_text = f"{1000.0 * delta:+.1f}"
                unit_text = r"\AA{} / m\AA"
            else:
                qc_text = f"{qc:.2f}"
                ensemble_text = f"{ensemble:.2f}"
                separate_text = f"{separate:.2f}"
                delta_text = f"{delta:+.2f}"
                unit_text = r"deg"
            handle.write(
                f"{class_name} & {unit_text} & {n_items:d} & {qc_text} & {ensemble_text} & "
                f"{separate_text} & {delta_text} \\\\\n"
            )
        handle.write(r"\bottomrule" + "\n")
        handle.write(r"\end{tabular}" + "\n")
    _write_summary(path.with_suffix(".summary.txt"), rows, title=title)


def _format_float(value: float, digits: int = 2) -> str:
    if not np.isfinite(value):
        return r"$\infty$"
    if value != 0.0 and (abs(value) >= 1.0e4 or abs(value) <= 1.0e-3):
        exponent = int(np.floor(np.log10(abs(value))))
        mantissa = value / (10.0**exponent)
        return rf"${mantissa:.2f}\times 10^{{{exponent}}}$"
    return f"{value:.{digits}f}"


def _write_glycine_typing_table(path: Path) -> None:
    rows = [
        ("Element and value windows", ROOT / "working/semiexp/glycine_ensemble/ensemble_manifest.json"),
        (r"Continuous synthon $Z_\mathrm{eff}$ types", ROOT / "working/semiexp/glycine_ensemble_synthon/ensemble_manifest.json"),
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write(r"\begin{tabular}{lrrrr}" + "\n")
        handle.write(r"\toprule" + "\n")
        handle.write(r"Class definition & Rank & Scaled condition & WRMS before & WRMS after \\" + "\n")
        handle.write(r"\midrule" + "\n")
        for label, manifest in rows:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            handle.write(
                f"{label} & {int(data['rank'])} & {_format_float(float(data['scaled_condition_number']), 2)} & "
                f"{float(data['weighted_rms_before']):.2f} & {float(data['weighted_rms_after']):.2f} \\\\\n"
            )
        handle.write(r"\bottomrule" + "\n")
        handle.write(r"\end{tabular}" + "\n")


def _write_glycine_threshold_table(path: Path) -> None:
    rows = _csv_rows(ROOT / "working/semiexp/glycine_synthon_threshold_scan_wide/synthon_threshold_scan.csv")
    with path.open("w", encoding="utf-8") as handle:
        handle.write(r"\begin{tabular}{rlrrrr}" + "\n")
        handle.write(r"\toprule" + "\n")
        handle.write(r"$Z_\mathrm{eff}$ threshold & Status & Rank & Cond. & WRMS & Min. matches \\" + "\n")
        handle.write(r"\midrule" + "\n")
        for row in rows:
            handle.write(
                f"{float(row['synthon_threshold']):.3f} & {row['acceptance_status']} & {int(row['rank'])} & "
                f"{_format_float(float(row['scaled_condition_number']), 2)} & "
                f"{float(row['weighted_rms_after']):.2f} & {int(row['min_matched_coordinates'])} \\\\\n"
            )
        handle.write(r"\bottomrule" + "\n")
        handle.write(r"\end{tabular}" + "\n")


def _write_glycine_classes_table(path: Path) -> None:
    rows = _csv_rows(ROOT / "working/semiexp/glycine_ensemble_synthon/ensemble_class_report.csv")
    with path.open("w", encoding="utf-8") as handle:
        handle.write(r"\begin{tabular}{llrrr}" + "\n")
        handle.write(r"\toprule" + "\n")
        handle.write(r"Class & Type & $Z_\mathrm{eff}$ prototype & Matches & Correction \\" + "\n")
        handle.write(r"\midrule" + "\n")
        for row in rows:
            name = row["class"].replace("_", "--")
            zeff = row.get("synthon_zeff", "").replace("|", ", ")
            handle.write(
                f"{name} & {row['kind']} & {zeff} & {int(float(row['matched_coordinates']))} & "
                f"{float(row['correction']):.5f} \\\\\n"
            )
        handle.write(r"\bottomrule" + "\n")
        handle.write(r"\end{tabular}" + "\n")


def _write_parent_only_fit_table(path: Path) -> None:
    data = json.loads((ROOT / "working/semiexp/anhydrides_parent_only/ensemble_manifest.json").read_text(encoding="utf-8"))
    with path.open("w", encoding="utf-8") as handle:
        handle.write(r"\begin{tabular}{lrrrr}" + "\n")
        handle.write(r"\toprule" + "\n")
        handle.write(r"Case & Parents & Obs. & Classes & Rank \\" + "\n")
        handle.write(r"\midrule" + "\n")
        for label in ("Maleic", "Phthalic", "Succinic"):
            handle.write(f"{label} parent & 1 & 2 & 1 & 1 \\\\\n")
        handle.write(
            f"Three-parent ensemble & 3 & 6 & {len(data['classes'])} & {int(data['rank'])} \\\\\n"
        )
        handle.write(r"\bottomrule" + "\n")
        handle.write(r"\end{tabular}" + "\n")


def _write_parent_only_classes_table(path: Path) -> None:
    rows = _csv_rows(ROOT / "working/semiexp/anhydrides_parent_only/ensemble_class_report.csv")
    with path.open("w", encoding="utf-8") as handle:
        handle.write(r"\begin{tabular}{lrrrr}" + "\n")
        handle.write(r"\toprule" + "\n")
        handle.write(r"Class & Matches & Molecules & Correction / \AA{} & $\sigma$ / \AA{} \\" + "\n")
        handle.write(r"\midrule" + "\n")
        for row in rows:
            class_name = row["class"].replace("_", r"\_")
            handle.write(
                f"{class_name} & {int(float(row['matched_coordinates']))} & "
                f"{int(float(row['molecule_count']))} & {float(row['correction']):+.5f} & "
                f"{float(row['sigma']):.5f} \\\\\n"
            )
        handle.write(r"\bottomrule" + "\n")
        handle.write(r"\end{tabular}" + "\n")


def _write_summary(path: Path, rows: list[Row], *, title: str) -> None:
    lines = [title]
    for molecule in sorted({row.molecule for row in rows}):
        items = [row for row in rows if row.molecule == molecule]
        stretch = [abs(row.delta) * 1000.0 for row in items if row.unit == "A"]
        bend = [abs(row.delta) for row in items if row.unit == "deg"]
        if stretch:
            lines.append(f"{molecule}: mean/max stretch ensemble-separate = {np.mean(stretch):.2f}/{np.max(stretch):.2f} mA")
        if bend:
            lines.append(f"{molecule}: mean/max bend ensemble-separate = {np.mean(bend):.3f}/{np.max(bend):.3f} deg")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    anhydrides = _geometry_rows(
        set_name="anhydrides",
        ensemble_job=ROOT / "examples/semiexp/anhydrides_parent_only/anhydrides_parent_only.mse-ensemble.toml",
        corrections_csv=ROOT / "working/semiexp/anhydrides_parent_only/ensemble_class_corrections.csv",
        separate_dirs={
            "maleic_anhydride": ROOT / "working/semiexp/maleic_anhydride",
            "phthalic_anhydride": ROOT / "working/semiexp/phthalic_anhydride",
            "succinic_anhydride": ROOT / "working/semiexp/succinic_anhydride_fixed_h",
        },
    )
    glycine = _geometry_rows(
        set_name="glycine",
        ensemble_job=ROOT / "examples/semiexp/glycine_ensemble/glycine_conformers_synthon.mse-ensemble.toml",
        corrections_csv=ROOT / "working/semiexp/glycine_ensemble_synthon/ensemble_class_corrections.csv",
        separate_dirs={
            "glycine_I": ROOT / "working/semiexp/glycine_I_cs_table5_repo_gic_fortran",
            "glycine_II": ROOT / "working/semiexp/glycine_II_dpcs3_table5_repo_gic_fortran",
        },
    )
    _write_detail_csv(ANALYSIS / "anhydrides_geometry_comparison.csv", anhydrides)
    _write_detail_csv(ANALYSIS / "glycine_geometry_comparison.csv", glycine)
    GENERATED.mkdir(parents=True, exist_ok=True)
    _write_table(
        GENERATED / "anhydrides_geometry_comparison.tex",
        anhydrides,
        title="Anhydride ensemble-vs-separate geometry comparison",
    )
    _write_table(
        GENERATED / "glycine_geometry_comparison.tex",
        glycine,
        title="Glycine ensemble-vs-separate geometry comparison",
    )
    _write_glycine_typing_table(GENERATED / "glycine_typing_comparison.tex")
    _write_glycine_threshold_table(GENERATED / "glycine_synthon_threshold_scan.tex")
    _write_glycine_classes_table(GENERATED / "glycine_synthon_classes.tex")
    _write_parent_only_fit_table(GENERATED / "anhydrides_parent_only_identifiability.tex")
    _write_parent_only_classes_table(GENERATED / "anhydrides_parent_only_classes.tex")


if __name__ == "__main__":
    main()
