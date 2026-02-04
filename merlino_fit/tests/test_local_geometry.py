import numpy as np

from survibfit.symmetry_local import match_local_geometry


def _unit(v):
    v = np.array(v, dtype=float)
    return v / np.linalg.norm(v)


def _match(label, vecs, tol=12.0):
    g = match_local_geometry(np.array([_unit(v) for v in vecs]), tol_deg=tol)
    return g.label


def test_geometry_5_tbp():
    vecs = [
        (0, 0, 1),
        (0, 0, -1),
        (1, 0, 0),
        (-0.5, np.sqrt(3) / 2.0, 0),
        (-0.5, -np.sqrt(3) / 2.0, 0),
    ]
    assert _match("tbp", vecs) == "tbp"


def test_geometry_6_octa():
    vecs = [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    ]
    assert _match("octa", vecs) == "octa"


def test_geometry_8_cube():
    vecs = [
        (1, 1, 1),
        (1, 1, -1),
        (1, -1, 1),
        (1, -1, -1),
        (-1, 1, 1),
        (-1, 1, -1),
        (-1, -1, 1),
        (-1, -1, -1),
    ]
    assert _match("cube", vecs) == "cube"


def test_geometry_8_dodeca():
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vecs = [
        (1, 0, phi),
        (-1, 0, phi),
        (1, 0, -phi),
        (-1, 0, -phi),
        (0, 1, phi),
        (0, -1, phi),
        (0, 1, -phi),
        (0, -1, -phi),
    ]
    assert _match("dodeca", vecs) == "dodeca"


def test_geometry_5_ambiguous():
    # Mix TBP and square-pyramidal to force ambiguity
    tbp = np.array(
        [
            (0, 0, 1),
            (0, 0, -1),
            (1, 0, 0),
            (-0.5, np.sqrt(3) / 2.0, 0),
            (-0.5, -np.sqrt(3) / 2.0, 0),
        ],
        dtype=float,
    )
    sp = np.array(
        [
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (0, 0, 1),
        ],
        dtype=float,
    )
    vecs = 0.5 * (tbp + sp)
    assert _match("ambiguous", vecs, tol=20.0) == "ambiguous"


def test_geometry_4_square_planar():
    vecs = [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
    ]
    assert _match("square_planar", vecs) == "square_planar"
