import numpy as np

from survibfit.transforms import cart_to_internal_grad, cart_to_internal_hess, cart_to_internal_coords, internal_to_cart_coords
from survibfit.pipeline import primitives_from_topology


def test_transform_stationary_consistency():
    # simple random geometry
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(4, 3))
    Z = np.array([6, 6, 6, 6])
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))

    ncart = coords.size
    Hx = np.eye(ncart)
    gx = np.zeros(ncart)

    Hs = cart_to_internal_hess(Hx, gx, coords, prims, include_curvature=False)
    # should be symmetric and PSD
    assert np.max(np.abs(Hs - Hs.T)) < 1e-8


def test_coord_roundtrip_small():
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(4, 3))
    Z = np.array([6, 6, 6, 6])
    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    q = cart_to_internal_coords(coords, prims, U=None)
    coords2 = internal_to_cart_coords(q, coords, prims, U=None, max_iter=5)
    assert np.max(np.abs(coords2 - coords)) < 1e-3
