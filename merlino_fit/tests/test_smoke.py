import numpy as np

from survibfit.pipeline import primitives_from_topology, eval_primitives
from survibfit.terms import generate_terms, eval_terms


def test_smoke_primitives_and_terms():
    # simple water-like geometry in au
    coords = np.array([
        [0.000000, 0.000000, 0.000000],
        [1.430000, 0.000000, 0.000000],
        [-0.358000, 1.330000, 0.000000],
    ])
    Z = np.array([8, 1, 1])

    prims = primitives_from_topology(coords, Z, linear_threshold=np.deg2rad(170.0))
    s = eval_primitives(prims, coords)
    assert s.size > 0

    nvib = min(3, s.size)
    q = s[:nvib]
    terms = generate_terms(nvib)
    basis_cfg = {i: {"mode": "poly", "params": {}} for i in range(nvib)}
    phi, dphi = eval_terms(q, terms, basis_cfg)
    assert phi.shape[0] == len(terms)
    assert dphi.shape[1] == len(terms)
