from __future__ import annotations

import csv
from pathlib import Path

from .contracts import IsotopologueObservation, RotationalConstants, VibrationalCorrection


CSV_FIELDS = (
    "label",
    "A_MHz",
    "B_MHz",
    "C_MHz",
    "delta_A_MHz",
    "delta_B_MHz",
    "delta_C_MHz",
    "correction_source",
    "substitutions",
)


def parse_substitutions(text: str) -> dict[int, int]:
    """Parse compact isotope substitutions like `2:13;5:18`."""
    result: dict[int, int] = {}
    text = (text or "").strip()
    if not text:
        return result
    for chunk in text.split(";"):
        atom_text, isotope_text = chunk.split(":", 1)
        atom_index = int(atom_text.strip())
        isotope_a = int(isotope_text.strip())
        if atom_index < 1:
            raise ValueError("Substitution atom indexes are one-based")
        result[atom_index] = isotope_a
    return result


def format_substitutions(substitutions: dict[int, int]) -> str:
    return ";".join(f"{atom}:{mass}" for atom, mass in sorted(substitutions.items()))


def read_observations_csv(path: Path) -> tuple[IsotopologueObservation, ...]:
    observations: list[IsotopologueObservation] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = set(CSV_FIELDS).difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing semiexp CSV columns: {', '.join(sorted(missing))}")
        for row in reader:
            observations.append(
                IsotopologueObservation(
                    label=str(row["label"]).strip(),
                    constants=RotationalConstants(
                        float(row["A_MHz"]),
                        float(row["B_MHz"]),
                        float(row["C_MHz"]),
                    ),
                    correction=VibrationalCorrection(
                        float(row["delta_A_MHz"] or 0.0),
                        float(row["delta_B_MHz"] or 0.0),
                        float(row["delta_C_MHz"] or 0.0),
                        source=str(row["correction_source"] or "unspecified"),
                    ),
                    substitutions=parse_substitutions(str(row["substitutions"] or "")),
                )
            )
    return tuple(observations)


def write_observations_csv(path: Path, observations: tuple[IsotopologueObservation, ...]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for obs in observations:
            writer.writerow(
                {
                    "label": obs.label,
                    "A_MHz": f"{obs.constants.A_MHz:.12g}",
                    "B_MHz": f"{obs.constants.B_MHz:.12g}",
                    "C_MHz": f"{obs.constants.C_MHz:.12g}",
                    "delta_A_MHz": f"{obs.correction.delta_A_MHz:.12g}",
                    "delta_B_MHz": f"{obs.correction.delta_B_MHz:.12g}",
                    "delta_C_MHz": f"{obs.correction.delta_C_MHz:.12g}",
                    "correction_source": obs.correction.source,
                    "substitutions": format_substitutions(obs.substitutions),
                }
            )
    return target


def corrected_constants_rows(
    observations: tuple[IsotopologueObservation, ...],
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for obs in observations:
        corrected = obs.corrected
        rows.append(
            {
                "label": obs.label,
                "A_e_MHz": corrected.A_MHz,
                "B_e_MHz": corrected.B_MHz,
                "C_e_MHz": corrected.C_MHz,
                "correction_source": obs.correction.source,
            }
        )
    return rows
