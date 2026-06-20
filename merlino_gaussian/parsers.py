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


def summarize_gaussian_log(path: Path) -> GaussianLogSummary:
    """Extract stable high-level markers from a Gaussian log/out file."""
    target = Path(path)
    text = target.read_text(encoding="utf-8", errors="replace")
    scf = tuple(float(match.group(1)) for match in SCF_RE.finditer(text))
    return GaussianLogSummary(
        path=target,
        normal_termination=NORMAL_TERMINATION in text,
        scf_energies_hartree=scf,
        standard_orientation_count=text.count(STANDARD_ORIENTATION),
        input_orientation_count=text.count(INPUT_ORIENTATION),
        scan_marker_count=text.lower().count("scan"),
        puckering_marker_count=text.count("QPck") + text.count("PhiP") + text.count("RPck"),
    )
