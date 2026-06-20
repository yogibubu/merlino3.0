from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


SCF_RE = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+]?\d+\.\d+)")
NORMAL_TERMINATION = "Normal termination of Gaussian"
STANDARD_ORIENTATION = "Standard orientation:"
INPUT_ORIENTATION = "Input orientation:"


@dataclass(frozen=True)
class GaussianLogSummary:
    path: Path
    normal_termination: bool
    scf_energies_hartree: tuple[float, ...]
    standard_orientation_count: int
    input_orientation_count: int
    scan_marker_count: int
    puckering_marker_count: int
    frequencies_cm: tuple[float, ...] = ()
    last_orientation: tuple[tuple[str, float, float, float], ...] = ()


def summarize_gaussian_log(path: Path) -> GaussianLogSummary:
    """Extract stable high-level markers from a Gaussian log/out file."""
    target = Path(path)
    text = target.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    scf = tuple(float(match.group(1)) for match in SCF_RE.finditer(text))
    return GaussianLogSummary(
        path=target,
        normal_termination=NORMAL_TERMINATION in text,
        scf_energies_hartree=scf,
        standard_orientation_count=text.count(STANDARD_ORIENTATION),
        input_orientation_count=text.count(INPUT_ORIENTATION),
        scan_marker_count=text.lower().count("scan"),
        puckering_marker_count=text.count("QPck") + text.count("PhiP") + text.count("RPck"),
        frequencies_cm=tuple(_parse_frequencies(text)),
        last_orientation=tuple(_parse_last_orientation(lines)),
    )


def _parse_frequencies(text: str) -> list[float]:
    values: list[float] = []
    for line in text.splitlines():
        if "Frequencies --" not in line:
            continue
        for token in line.split("--", 1)[1].split():
            try:
                values.append(float(token))
            except ValueError:
                continue
    return values


def _parse_last_orientation(lines: list[str]) -> list[tuple[str, float, float, float]]:
    start = -1
    for idx, line in enumerate(lines):
        if STANDARD_ORIENTATION in line or INPUT_ORIENTATION in line:
            start = idx
    if start < 0:
        return []
    atoms: list[tuple[str, float, float, float]] = []
    dash_count = 0
    for line in lines[start + 1:]:
        if set(line.strip()) == {"-"}:
            dash_count += 1
            if dash_count >= 3:
                break
            continue
        if dash_count < 2:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            atomic_number = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
        except ValueError:
            continue
        atoms.append((str(atomic_number), x, y, z))
    return atoms
