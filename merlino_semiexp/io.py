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
    "sigma_A_MHz",
    "sigma_B_MHz",
    "sigma_C_MHz",
)

REQUIRED_CSV_FIELDS = (
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
        missing = set(REQUIRED_CSV_FIELDS).difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing semiexp CSV columns: {', '.join(sorted(missing))}")
        for row in reader:
            weights = _weights_from_sigmas(row)
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
                    weights=weights,
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
                    "sigma_A_MHz": _sigma_text(obs.weights.A_MHz) if obs.weights else "",
                    "sigma_B_MHz": _sigma_text(obs.weights.B_MHz) if obs.weights else "",
                    "sigma_C_MHz": _sigma_text(obs.weights.C_MHz) if obs.weights else "",
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


def _weights_from_sigmas(row: dict[str, str]) -> RotationalConstants | None:
    keys = ("sigma_A_MHz", "sigma_B_MHz", "sigma_C_MHz")
    raw = tuple(str(row.get(key, "") or "").strip() for key in keys)
    if not any(raw):
        return None
    if not all(raw):
        raise ValueError("Semiexp sigma columns must be all present or all empty")
    sigmas = tuple(float(item) for item in raw)
    if any(sigma <= 0.0 for sigma in sigmas):
        raise ValueError("Semiexp sigma columns must be positive when provided")
    return RotationalConstants(*(1.0 / (sigma * sigma) for sigma in sigmas))


def _sigma_text(weight: float) -> str:
    if weight <= 0.0:
        return ""
    return f"{(1.0 / weight) ** 0.5:.12g}"
