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
class QMParameterPredicate:
    """Weighted QM prior/predicate for a generated GIC parameter."""

    label_pattern: str
    value: float
    sigma: float
    source: str = "qm"

    @property
    def weight(self) -> float:
        if self.sigma <= 0.0:
            raise ValueError("QM predicate sigma must be positive")
        return 1.0 / (self.sigma * self.sigma)


@dataclass(frozen=True)
class SemiexperimentalFitRequest:
    initial_geometry: Path
    observations: tuple[IsotopologueObservation, ...]
    fixed_parameters: tuple[str, ...] = ()
    observable: str = "moments"
    rotational_components: str = "auto"
    qm_predicates: tuple[QMParameterPredicate, ...] = ()

    def validate(self) -> None:
        if not self.observations:
            raise ValueError("Semiexperimental fit needs at least one isotopologue")
        labels = [item.label for item in self.observations]
        if len(set(labels)) != len(labels):
            raise ValueError("Duplicate isotopologue labels are not allowed")
        if self.observable not in {"moments", "rotational_constants", "auto"}:
            raise ValueError("observable must be moments, rotational_constants or auto")
        if self.rotational_components not in {"auto", "ABC", "AB", "AC", "BC"}:
            raise ValueError("rotational_components must be auto, ABC, AB, AC or BC")
        for predicate in self.qm_predicates:
            if not predicate.label_pattern.strip():
                raise ValueError("QM predicate label pattern cannot be empty")
            if predicate.sigma <= 0.0:
                raise ValueError("QM predicate sigma must be positive")
