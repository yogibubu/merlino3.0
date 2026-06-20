from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class ForceFieldSource:
    """Normalized source descriptor for a Gaussian quartic force field."""

    path: Path
    source_type: str = "gaussian-log"
    checksum: str | None = None


@dataclass(frozen=True)
class DavidsonSettings:
    """Numerical settings for future large VCI diagonalization."""

    n_roots: int = 10
    max_subspace: int = 80
    max_iter: int = 200
    convergence: float = 1.0e-8
    preconditioner: str = "diagonal"

    def validate(self) -> None:
        if self.n_roots < 1:
            raise ValueError("Davidson needs at least one root")
        if self.max_subspace < self.n_roots:
            raise ValueError("Davidson max_subspace must be >= n_roots")
        if self.max_iter < 1:
            raise ValueError("Davidson max_iter must be positive")
        if self.convergence <= 0.0:
            raise ValueError("Davidson convergence must be positive")


@dataclass(frozen=True)
class VCIRequest:
    """High-level VPT2/VCI request before backend-specific input is written."""

    force_field: ForceFieldSource
    max_quanta: int = 6
    basis_energy_cutoff_cm: float | None = None
    davidson: DavidsonSettings = field(default_factory=DavidsonSettings)

    def validate(self) -> None:
        if self.max_quanta < 1:
            raise ValueError("VCI max_quanta must be positive")
        if self.basis_energy_cutoff_cm is not None and self.basis_energy_cutoff_cm <= 0.0:
            raise ValueError("VCI basis energy cutoff must be positive")
        self.davidson.validate()


@dataclass(frozen=True)
class VPT2VCIInventory:
    """Current Fortran VPT2/VCI code discovered in the Merlino tree."""

    legacy_vci1d: Path | None
    notes: tuple[str, ...]


def inventory_legacy_fortran(repo_root: Path) -> VPT2VCIInventory:
    """Record what VPT2/VCI Fortran exists before implementing new kernels."""
    root = Path(repo_root)
    legacy_vci1d = root / "fortran" / "legacy" / "vibrational" / "vci1d.f"
    notes = [
        "Current tree contains a legacy one-dimensional VCI driver.",
        "Full quartic-field VPT2/VCI backend and Davidson diagonalizer are not active yet.",
    ]
    return VPT2VCIInventory(
        legacy_vci1d=legacy_vci1d if legacy_vci1d.exists() else None,
        notes=tuple(notes),
    )
