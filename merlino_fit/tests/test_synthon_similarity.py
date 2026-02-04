from pathlib import Path

from survibfit.synthon_similarity import compare_molecules, compare_against_library


def test_similarity_is_high_for_identical_molecule():
    xyz = Path(__file__).resolve().parent / "data" / "ch4.xyz"
    out = compare_molecules(xyz, xyz, covariance_mode="full")
    assert out["similarity_exp_minus_db"] > 0.999
    assert out["bhattacharyya_distance"] < 1.0e-6


def test_similarity_distinguishes_different_molecules():
    root = Path(__file__).resolve().parent / "data"
    ch4 = root / "ch4.xyz"
    co2 = root / "co2.xyz"
    same = compare_molecules(ch4, ch4, covariance_mode="full")
    diff = compare_molecules(ch4, co2, covariance_mode="full")
    assert diff["similarity_exp_minus_db"] < same["similarity_exp_minus_db"]


def test_library_ranking_is_sorted():
    root = Path(__file__).resolve().parent / "data"
    ch4 = root / "ch4.xyz"
    co2 = root / "co2.xyz"
    c4 = root / "c4_chain.xyz"
    out = compare_against_library(ch4, [ch4, co2, c4], covariance_mode="full")
    assert out["library_size"] == 2
    sims = [row["similarity_combined"] for row in out["ranking"]]
    assert sims == sorted(sims, reverse=True)


def test_ring_comparison_penalizes_ring_mismatch():
    root = Path(__file__).resolve().parent / "data"
    c5h5 = root / "c5h5.xyz"
    ch4 = root / "ch4.xyz"

    same = compare_molecules(c5h5, c5h5, covariance_mode="full")
    diff = compare_molecules(c5h5, ch4, covariance_mode="full")

    assert same["ring_similarity_exp_minus_db"] > 0.999
    assert diff["ring_similarity_exp_minus_db"] < same["ring_similarity_exp_minus_db"]
    assert diff["similarity_combined"] < same["similarity_combined"]
