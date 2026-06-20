from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement, product

import numpy as np


@dataclass(frozen=True)
class QuarticForceField:
    """Dimensionless normal-coordinate force field used by the VCI solver."""

    harmonic_frequencies_cm: np.ndarray
    cubic_cm: dict[tuple[int, int, int], float]
    quartic_cm: dict[tuple[int, int, int, int], float]


@dataclass(frozen=True)
class VCIResult:
    basis: tuple[tuple[int, ...], ...]
    energies_cm: np.ndarray
    eigenvectors: np.ndarray


def generate_vibrational_basis(n_modes: int, max_quanta: int) -> tuple[tuple[int, ...], ...]:
    """Generate product harmonic-oscillator states with total quanta cutoff."""
    if n_modes < 1 or max_quanta < 0:
        raise ValueError("n_modes must be positive and max_quanta non-negative")
    states = [state for state in product(range(max_quanta + 1), repeat=n_modes) if sum(state) <= max_quanta]
    return tuple(sorted(states, key=lambda s: (sum(s), s)))


def _x_matrix_power(max_n: int, power: int) -> np.ndarray:
    x = np.zeros((max_n + 1, max_n + 1), dtype=float)
    for n in range(max_n + 1):
        if n + 1 <= max_n:
            x[n + 1, n] = np.sqrt(n + 1) / np.sqrt(2.0)
        if n - 1 >= 0:
            x[n - 1, n] = np.sqrt(n) / np.sqrt(2.0)
    out = np.eye(max_n + 1)
    for _ in range(power):
        out = out @ x
    return out


def _mode_powers(indices: tuple[int, ...], n_modes: int) -> tuple[int, ...]:
    powers = [0] * n_modes
    for idx in indices:
        if idx < 0 or idx >= n_modes:
            raise ValueError("Force-field mode index out of range")
        powers[idx] += 1
    return tuple(powers)


def _term_element(left: tuple[int, ...], right: tuple[int, ...], xpowers: tuple[np.ndarray, ...]) -> float:
    value = 1.0
    for mode, power_matrix in enumerate(xpowers):
        value *= power_matrix[left[mode], right[mode]]
        if value == 0.0:
            break
    return value


def build_vci_hamiltonian(force_field: QuarticForceField, max_quanta: int) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Build a small dense VCI Hamiltonian in cm^-1.

    Large spaces will later use the same matrix-element code behind Davidson.
    """
    freqs = np.asarray(force_field.harmonic_frequencies_cm, dtype=float)
    if np.any(freqs <= 0.0):
        raise ValueError("Harmonic frequencies must be positive")
    n_modes = len(freqs)
    basis = generate_vibrational_basis(n_modes, max_quanta)
    max_n = max(max(state) for state in basis)
    operator_max_n = max_n + 4
    powers = {p: _x_matrix_power(operator_max_n, p) for p in range(5)}
    h = np.zeros((len(basis), len(basis)), dtype=float)

    for i, state in enumerate(basis):
        h[i, i] = float(np.dot(freqs, np.array(state, dtype=float) + 0.5))

    terms: list[tuple[tuple[int, ...], float]] = []
    terms.extend((tuple(sorted(k)), v) for k, v in force_field.cubic_cm.items())
    terms.extend((tuple(sorted(k)), v) for k, v in force_field.quartic_cm.items())
    for indices, coeff in terms:
        mode_powers = _mode_powers(indices, n_modes)
        xpowers = tuple(powers[p] for p in mode_powers)
        for row, left in enumerate(basis):
            for col, right in enumerate(basis[: row + 1]):
                element = coeff * _term_element(left, right, xpowers)
                if element:
                    h[row, col] += element
                    if row != col:
                        h[col, row] += element
    return h, basis


def solve_vci(force_field: QuarticForceField, max_quanta: int, n_roots: int | None = None) -> VCIResult:
    hamiltonian, basis = build_vci_hamiltonian(force_field, max_quanta)
    eig, vec = np.linalg.eigh((hamiltonian + hamiltonian.T) * 0.5)
    if n_roots is not None:
        eig = eig[:n_roots]
        vec = vec[:, :n_roots]
    return VCIResult(basis=basis, energies_cm=eig, eigenvectors=vec)


def zero_anharmonic_force_field(frequencies_cm: np.ndarray) -> QuarticForceField:
    return QuarticForceField(np.asarray(frequencies_cm, dtype=float), {}, {})


def empty_symmetric_terms(n_modes: int, order: int) -> dict[tuple[int, ...], float]:
    """Return all same-index combinations initialized to zero for tests/UI forms."""
    if order not in {3, 4}:
        raise ValueError("Only cubic and quartic terms are supported")
    return {tuple(idx): 0.0 for idx in combinations_with_replacement(range(n_modes), order)}
