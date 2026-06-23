import numpy as np

from survibfit.symmetry_global import _group_label, irrep_characters_for_operations


def test_group_cnh():
    labels = ["C3z^1", "sigma_xy", "C2_xy_3_0"]
    pg = _group_label([(lab, None, 0.0) for lab in labels])
    assert pg.startswith("C3h")


def test_group_dnd():
    labels = ["C4z^1", "C2_xy_4_0", "S4"]
    pg = _group_label([(lab, None, 0.0) for lab in labels])
    assert pg.startswith("D4d")


def test_group_sn():
    labels = ["S6"]
    pg = _group_label([(lab, None, 0.0) for lab in labels])
    assert pg.startswith("S6")


def test_group_c3h():
    labels = ["C3z^1", "sigma_xy"]
    pg = _group_label([(lab, None, 0.0) for lab in labels])
    assert pg.startswith("C3h")


def test_group_d3d():
    labels = ["C3z^1", "C2_xy_3_0", "S6"]
    pg = _group_label([(lab, None, 0.0) for lab in labels])
    assert pg.startswith("D3d")


def test_irrep_character_tables_cover_common_point_group_families():
    cases = [
        ("C3", ["E", "C3z^1", "C3z^2"], {"A", "E1"}),
        ("C3v", ["E", "C3z^1", "C3z^2", "sigma_v_3_0"], {"A1", "A2", "E1"}),
        ("C3h", ["E", "C3z^1", "C3z^2", "sigma_xy"], {"Ag", "Au", "E1g", "E1u"}),
        ("D3", ["E", "C3z^1", "C3z^2", "C2_xy_3_0"], {"A1", "A2", "E1"}),
        ("D3d", ["E", "C3z^1", "C3z^2", "C2_xy_3_0", "i", "S6"], {"A1g", "A1u", "A2g", "A2u", "E1g", "E1u"}),
        ("D4d", ["E", "C4z^1", "C4z^3", "C2_xy_4_0", "S4"], {"A1'", "A1''", "A2'", "A2''", "B1'", "B1''", "B2'", "B2''", "E1'", "E1''"}),
        ("Td", ["E", "C3_t", "C2_t", "S4", "sigma_xy"], {"A1", "A2", "E", "T1", "T2"}),
        ("Oh", ["E", "C4_o", "C3_o", "C2_o", "i", "sigma_xy"], {"A1g", "A1u", "A2g", "A2u", "Eg", "Eu", "T1g", "T1u", "T2g", "T2u"}),
        ("Ih", ["E", "C5_i", "C5_i2", "C3_i", "C2_i", "i"], {"Ag", "Au", "T1g", "T1u", "T2g", "T2u", "Gg", "Gu", "Hg", "Hu"}),
    ]
    for point_group, labels, expected in cases:
        irreps = irrep_characters_for_operations(labels, point_group=point_group)
        names = {name for name, chars in irreps}
        assert expected <= names
        assert all(chars.shape == (len(labels),) for _name, chars in irreps)


def test_group_cinfv():
    labels = ["C3z^1", "C5z^1", "sigma_v_6_0", "sigma_v_6_1", "sigma_v_6_2"]
    pg = _group_label([(lab, None, 0.0) for lab in labels])
    assert pg in ("Cinfv", "Dinfh")


def test_op_filter_excludes_mirrors():
    labels = ["C3z^1", "sigma_xy", "sigma_v_6_0"]
    from survibfit.symmetry_global import symmetry_elements_from_geometry, _group_label, orient_coords
    # build dummy coords for two atoms to trigger operations list filtering
    coords = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]])
    symbols = ["H", "H"]
    coords_oriented = orient_coords(coords)
    elements, classes, perms = symmetry_elements_from_geometry(
        symbols,
        coords_oriented,
        tol=1e-3,
        max_n=6,
        op_filter=lambda lab: not lab.startswith("sigma"),
    )
    pg = _group_label(elements)
    assert "sigma" not in " ".join([e[0] for e in elements])


def test_no_radial_filter_does_not_crash():
    from survibfit.symmetry_global import symmetry_elements_from_geometry, orient_coords
    coords = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]])
    symbols = ["H", "H"]
    coords_oriented = orient_coords(coords)
    elements, _, _ = symmetry_elements_from_geometry(
        symbols,
        coords_oriented,
        tol=1e-3,
        max_n=6,
        enforce_radial_filter=False,
    )
    assert elements


def test_co2_linear_dinfh():
    from pathlib import Path
    from survibfit.modify_geom import read_xyz
    from topology.elements import atomic_number
    from survibfit.symmetry_global import orient_coords, symmetry_elements_from_geometry, _group_label

    xyz = Path(__file__).resolve().parent / "data" / "co2.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    symbols = [a for a in atoms]
    coords_oriented = orient_coords(coords_ang, weights=Z)
    elements, classes, permutations = symmetry_elements_from_geometry(
        symbols, coords_oriented, tol=1.0e-3, max_n=10
    )
    pg = _group_label(elements)
    assert pg in ("Dinfh", "Cinfv", "D10h")


def test_benzene_not_overclassified_as_d8h():
    from survibfit.symmetry_global import orient_coords, symmetry_elements_from_geometry, _group_label

    coords = np.array(
        [
            [1.396, 0.000, 0.000],
            [0.698, 1.209, 0.000],
            [-0.698, 1.209, 0.000],
            [-1.396, 0.000, 0.000],
            [-0.698, -1.209, 0.000],
            [0.698, -1.209, 0.000],
            [2.479, 0.000, 0.000],
            [1.240, 2.148, 0.000],
            [-1.240, 2.148, 0.000],
            [-2.479, 0.000, 0.000],
            [-1.240, -2.148, 0.000],
            [1.240, -2.148, 0.000],
        ],
        dtype=float,
    )
    symbols = ["C"] * 6 + ["H"] * 6
    coords_oriented = orient_coords(coords)
    elements, _, _ = symmetry_elements_from_geometry(
        symbols,
        coords_oriented,
        tol=1.0e-3,
        max_n=8,
        auto_max_n=True,
    )
    pg = _group_label(elements)
    assert pg == "D6h"


def test_pyridine_is_c2v():
    from pathlib import Path
    from survibfit.modify_geom import read_xyz
    from topology.elements import atomic_number
    from survibfit.symmetry_global import orient_coords, symmetry_elements_from_geometry, _group_label

    xyz = Path(__file__).resolve().parents[1] / "pyridine_in.xyz"
    atoms, coords_ang, _ = read_xyz(xyz)
    Z = np.array([atomic_number(a) for a in atoms], dtype=int)
    symbols = [a for a in atoms]
    coords_oriented = orient_coords(coords_ang, weights=Z)
    elements, _, _ = symmetry_elements_from_geometry(
        symbols, coords_oriented, tol=1.0e-3, max_n=10
    )
    pg = _group_label(elements)
    assert pg == "C2v"
