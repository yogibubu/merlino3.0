from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HessianInput:
    """Canonical Cartesian harmonic input for Merlino GF.

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
