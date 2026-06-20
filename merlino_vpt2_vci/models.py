from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HessianInput:
    """Merlino canonical Cartesian harmonic input.

    Coordinates are in bohr, masses in amu, and the Cartesian Hessian is a full
    symmetric matrix in Eh/bohr^2. No file format is implied.
    """

    atomic_numbers: np.ndarray
    cartesian_coordinates_bohr: np.ndarray
    masses_amu: np.ndarray
    cartesian_hessian: np.ndarray
    harmonic_frequencies_cm: np.ndarray
    source: str = "merlino"

    def validate(self) -> None:
        natoms = len(self.atomic_numbers)
        if self.cartesian_coordinates_bohr.shape != (natoms, 3):
            raise ValueError("Cartesian coordinates must have shape (natoms, 3)")
        if self.masses_amu.shape != (natoms,):
            raise ValueError("Masses must have shape (natoms,)")
        expected = (3 * natoms, 3 * natoms)
        if self.cartesian_hessian.shape != expected:
            raise ValueError(f"Cartesian Hessian must have shape {expected}")
        if not np.allclose(self.cartesian_hessian, self.cartesian_hessian.T):
            raise ValueError("Cartesian Hessian must be symmetric")


@dataclass(frozen=True)
class AnharmonicInput:
    """Merlino canonical normal-coordinate anharmonic input."""

    harmonic_frequencies_cm: np.ndarray
    anharmonic_frequencies_cm: np.ndarray
    cubic_cm: dict[tuple[int, int, int], float]
    quartic_cm: dict[tuple[int, int, int, int], float]
    source: str = "merlino"

    def validate(self) -> None:
        if self.harmonic_frequencies_cm.ndim != 1:
            raise ValueError("Harmonic frequencies must be a one-dimensional array")
        if self.anharmonic_frequencies_cm.ndim != 1:
            raise ValueError("Anharmonic frequencies must be a one-dimensional array")
