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
    """Current VPT2/VCI backend status discovered in the Merlino4 tree."""

    harmonic_internal_source: Path | None
    gdv_vci_driver_source: Path | None
    gdv_davidson_source: Path | None
    active_fortran_sources: tuple[Path, ...]
    davidson_backend: Path | None
    notes: tuple[str, ...]


def inventory_vpt2_vci_backends(repo_root: Path) -> VPT2VCIInventory:
    """Record active VPT2/VCI kernels available in Merlino4."""
    root = Path(repo_root)
    harmonic_internal = root / "fortran" / "harmonic_internal" / "gf.f"
    gdv_root = Path.home() / "gdv_j32p" / "gdv"
    gdv_vci_driver = gdv_root / "l717.F"
    gdv_davidson = gdv_root / "utilnz.F"
    source_dir = root / "fortran" / "vpt2_vci"
    sources = tuple(sorted(source_dir.glob("*.f"))) if source_dir.exists() else ()
    davidson = source_dir / "davidson_core.f"
    notes = [
        "Harmonic internal-coordinate GF analysis is available through gf.f.",
        "GDV l717.F contains the current VPT2/VCI driver decks.",
        "GDV utilnz.F is retained only as historical context; Merlino4 Davidson is independent.",
        "Merlino4 has independent Python and Fortran77 GF/VCI/Davidson cores; Gaussian QFF tensor promotion is still being expanded.",
    ]
    return VPT2VCIInventory(
        harmonic_internal_source=harmonic_internal if harmonic_internal.exists() else None,
        gdv_vci_driver_source=gdv_vci_driver if gdv_vci_driver.exists() else None,
        gdv_davidson_source=gdv_davidson if gdv_davidson.exists() else None,
        active_fortran_sources=sources,
        davidson_backend=davidson if davidson.exists() else None,
        notes=tuple(notes),
    )
