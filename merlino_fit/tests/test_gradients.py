import numpy as np

from survibfit.bmat import dihedral_grad, oop_grad, linear_grad
from survibfit.geometry import dihedral, oop, linear_components
from survibfit.bmat import finite_diff_grad


def _rand_coords(n=4):
    rng = np.random.default_rng(123)
    return rng.normal(size=(n, 3))


def test_dihedral_grad_matches_fd():
    coords = _rand_coords(4)
    g_ana = dihedral_grad(0, 1, 2, 3, coords)
    g_fd = finite_diff_grad(lambda c: dihedral(0, 1, 2, 3, c), coords, h=1e-6)
    assert np.max(np.abs(g_ana - g_fd)) < 1e-4


def test_oop_grad_matches_fd():
    coords = _rand_coords(4)
    g_ana = oop_grad(0, 1, 2, 3, coords)
    g_fd = finite_diff_grad(lambda c: oop(0, 1, 2, 3, c), coords, h=1e-6)
    assert np.max(np.abs(g_ana - g_fd)) < 1e-4


def test_linear_grad_matches_fd():
    coords = _rand_coords(3)
    g_ana = linear_grad(0, 1, 2, coords, mode=-1)
    g_fd = finite_diff_grad(lambda c: linear_components(0, 1, 2, c)[0], coords, h=1e-6)
    assert np.max(np.abs(g_ana - g_fd)) < 1e-4
