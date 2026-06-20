from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class RotationalConstants:
    """Rotational constants in MHz."""

    A_MHz: float
    B_MHz: float
    C_MHz: float

    def as_tuple(self) -> tuple[float, float, float]:
        return (self.A_MHz, self.B_MHz, self.C_MHz)


@dataclass(frozen=True)
class VibrationalCorrection:
    """QM vibrational correction Delta B_vib in MHz."""

    delta_A_MHz: float = 0.0
    delta_B_MHz: float = 0.0
    delta_C_MHz: float = 0.0
    source: str = "unspecified"

    def as_tuple(self) -> tuple[float, float, float]:
        return (self.delta_A_MHz, self.delta_B_MHz, self.delta_C_MHz)


@dataclass(frozen=True)
class CorrectedRotationalConstants:
    """Experimental constants corrected to semiexperimental equilibrium values."""

    observed: RotationalConstants
    correction: VibrationalCorrection

    @property
    def equilibrium(self) -> RotationalConstants:
        a, b, c = self.observed.as_tuple()
        da, db, dc = self.correction.as_tuple()
        return RotationalConstants(a - da, b - db, c - dc)


@dataclass(frozen=True)
class IsotopologueObservation:
    label: str
    constants: RotationalConstants
    substitutions: dict[int, int] = field(default_factory=dict)
    correction: VibrationalCorrection = field(default_factory=VibrationalCorrection)
    weights: RotationalConstants | None = None

    @property
    def corrected(self) -> RotationalConstants:
        return CorrectedRotationalConstants(self.constants, self.correction).equilibrium


@dataclass(frozen=True)
class SemiexperimentalFitRequest:
    initial_geometry: Path
    observations: tuple[IsotopologueObservation, ...]
    fixed_parameters: tuple[str, ...] = ()

    def validate(self) -> None:
        if not self.observations:
            raise ValueError("Semiexperimental fit needs at least one isotopologue")
        labels = [item.label for item in self.observations]
        if len(set(labels)) != len(labels):
            raise ValueError("Duplicate isotopologue labels are not allowed")
