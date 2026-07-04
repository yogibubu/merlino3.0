from __future__ import annotations

import csv
import json
from pathlib import Path

from merlino_semiexp import build_reference_assisted_geometry, load_se_reference_library, search_reference_library


def test_reference_library_search_ranks_closest_geometry_first(tmp_path):
    query = tmp_path / "query_water.xyz"
    query.write_text(
        "\n".join(
            [
                "3",
                "query water",
                "O 0.000000 0.000000 0.000000",
                "H 0.958400 0.000000 0.000000",
                "H -0.239000 0.927000 0.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    library = tmp_path / "library"
    xyz_dir = library / "xyz"
    xyz_dir.mkdir(parents=True)
    (xyz_dir / "water_ref.xyz").write_text(query.read_text(encoding="utf-8"), encoding="utf-8")
    (xyz_dir / "methane_ref.xyz").write_text(
        "\n".join(
            [
                "5",
                "methane reference",
                "C 0.000000 0.000000 0.000000",
                "H 0.629000 0.629000 0.629000",
                "H -0.629000 -0.629000 0.629000",
                "H -0.629000 0.629000 -0.629000",
                "H 0.629000 -0.629000 -0.629000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with (library / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["slug", "name", "atoms", "level", "path"])
        writer.writerow(["water_ref", "Water reference", 3, "SE", "xyz/water_ref.xyz"])
        writer.writerow(["methane_ref", "Methane reference", 5, "SE", "xyz/methane_ref.xyz"])

    result = search_reference_library(query, library_root=library, top_k=2, outdir=tmp_path / "out")

    assert result.matches[0].slug == "water_ref"
    assert result.matches[0].similarity_combined > result.matches[1].similarity_combined
    assert (tmp_path / "out" / "reference_matches.csv").exists()
    payload = json.loads((tmp_path / "out" / "reference_matches.json").read_text(encoding="utf-8"))
    assert payload["matches"][0]["slug"] == "water_ref"
    assert payload["settings"]["library_size"] == 2


def test_reference_assisted_geometry_transfers_supported_fragments(tmp_path):
    query, library = _write_tiny_reference_library(tmp_path)

    result = build_reference_assisted_geometry(
        query,
        library_root=library,
        top_library_matches=2,
        max_fragment_matches=2,
        min_fragment_support=1,
        zeff_threshold=0.2,
        apply_kinds=("bond", "angle"),
        outdir=tmp_path / "assisted",
    )

    assert result.targets
    assert result.rms_target_residual_final <= result.rms_target_residual_initial
    assert (tmp_path / "assisted" / "reference_assisted_geometry.xyz").exists()
    assert (tmp_path / "assisted" / "multiclasses_fragment_summary.csv").exists()


def test_default_se_reference_library_manifest_is_available():
    references = load_se_reference_library()

    assert len(references) >= 100
    assert all(ref.path.is_file() for ref in references)
    assert any("maleic" in ref.slug for ref in references)


def _write_tiny_reference_library(tmp_path: Path) -> tuple[Path, Path]:
    query = tmp_path / "query_water.xyz"
    query.write_text(
        "\n".join(
            [
                "3",
                "query water",
                "O 0.000000 0.000000 0.000000",
                "H 0.958400 0.000000 0.000000",
                "H -0.239000 0.927000 0.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    library = tmp_path / "library"
    xyz_dir = library / "xyz"
    xyz_dir.mkdir(parents=True)
    (xyz_dir / "water_ref.xyz").write_text(
        "\n".join(
            [
                "3",
                "water reference",
                "O 0.000000 0.000000 0.000000",
                "H 0.970000 0.000000 0.000000",
                "H -0.241890 0.938224 0.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (xyz_dir / "methane_ref.xyz").write_text(
        "\n".join(
            [
                "5",
                "methane reference",
                "C 0.000000 0.000000 0.000000",
                "H 0.629000 0.629000 0.629000",
                "H -0.629000 -0.629000 0.629000",
                "H -0.629000 0.629000 -0.629000",
                "H 0.629000 -0.629000 -0.629000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with (library / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["slug", "name", "atoms", "level", "path"])
        writer.writerow(["water_ref", "Water reference", 3, "SE", "xyz/water_ref.xyz"])
        writer.writerow(["methane_ref", "Methane reference", 5, "SE", "xyz/methane_ref.xyz"])
    return query, library
