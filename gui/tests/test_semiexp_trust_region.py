import numpy as np
import pytest

from merlino_semiexp.fit import (
    _accepted_trust_update,
    _predicted_reduction,
    _rejected_trust_update,
    _svd_trust_region_lm_step,
)


def test_svd_trust_region_returns_gauss_newton_inside_radius():
    jac = np.array([[2.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=float)
    residual = np.array([2.0, 1.0, 2.0], dtype=float)

    result = _svd_trust_region_lm_step(jac, residual, damping=1.0e-8, trust_radius=10.0)
    expected = np.linalg.lstsq(jac, residual, rcond=None)[0]

    assert result.solver == "svd_more_hebden_trust_region"
    assert result.on_boundary is False
    assert result.step == pytest.approx(expected)


def test_svd_trust_region_solves_boundary_step_when_radius_is_active():
    jac = np.array([[10.0, 0.0], [0.0, 1.0]], dtype=float)
    residual = np.array([100.0, 1.0], dtype=float)
    radius = 0.5

    result = _svd_trust_region_lm_step(jac, residual, damping=1.0e-8, trust_radius=radius)
    reduction = _predicted_reduction(
        residual,
        jac,
        result.step,
        scale=1.0,
        current_objective=0.5 * float(residual @ residual),
    )

    assert result.on_boundary is True
    assert result.shift > 0.0
    assert np.linalg.norm(result.step) == pytest.approx(radius, rel=1.0e-8, abs=1.0e-10)
    assert reduction > 0.0


def test_svd_trust_region_is_rank_revealing_for_deficient_jacobian():
    jac = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]], dtype=float)
    residual = np.array([1.0, 2.0, 3.0], dtype=float)

    result = _svd_trust_region_lm_step(jac, residual, damping=0.0, trust_radius=10.0)

    assert np.all(np.isfinite(result.step))
    assert result.step == pytest.approx([0.5, 0.5])
    assert result.on_boundary is False


def test_trust_radius_updates_expand_only_for_good_boundary_steps():
    damping, radius = _accepted_trust_update(
        damping=1.0e-4,
        trust_radius=0.1,
        ratio=0.9,
        scale=1.0,
        step_norm=0.095,
        max_step=1.0,
    )

    assert damping < 1.0e-4
    assert radius > 0.1


def test_trust_radius_updates_contract_after_rejection():
    damping, radius = _rejected_trust_update(damping=0.0, trust_radius=0.1, max_step=1.0)

    assert damping > 0.0
    assert radius < 0.1
