from __future__ import annotations

import json
import hashlib
from pathlib import Path

from merlino_gic import run_gicforge
from merlino_gic.gic_symmetry import write_gic_symmetry_files


def test_run_gicforge_collects_outputs_and_manifest(tmp_path):
    executable = tmp_path / "fake_gicforge.sh"
    executable.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -e",
                "printf 'legacy report\\n' > provout",
                "printf 'gaussian input\\n' > gauin",
                "printf 'vpt2 input\\n' > VPT2in",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    (tmp_path / "provin").write_text("input\n", encoding="utf-8")

    result = run_gicforge(tmp_path, executable=executable)

    assert result.files["provout"] == tmp_path / "provout"
    assert result.files["gicforge.out"].read_text(encoding="utf-8") == "legacy report\n"
    assert result.files["gauin"].read_text(encoding="utf-8") == "gaussian input\n"
    manifest = json.loads(result.manifest.read_text(encoding="utf-8"))
    assert manifest["workflow"] == "gicforge"
    assert manifest["outputs"]["gauin"] == str(tmp_path / "gauin")
    assert "provin" in manifest["inputs"]


def test_gic_symmetry_postcheck_is_byte_deterministic(tmp_path):
    (tmp_path / "xyzin").write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.757000 0.000000 0.586000",
                "H -0.757000 0.000000 0.586000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "gauin").write_text(
        "\n".join(
            [
                "%chk=water.chk",
                "",
                "0 1",
                "",
                " S1=[ 0.70710678*R(  1,  2)+0.70710678*R(  1,  3)]",
                " S2=[ 0.70710678*R(  1,  2)-0.70710678*R(  1,  3)]",
                " A1=[ 1.00000000*A(  2,  1,  3)]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    outputs = ("gauin.symm", "gicsym", "gic_symmetry_diagnostics.json")
    write_gic_symmetry_files(tmp_path)
    first = {name: hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() for name in outputs}
    write_gic_symmetry_files(tmp_path)
    second = {name: hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() for name in outputs}

    assert first == second
    diagnostics = json.loads((tmp_path / "gic_symmetry_diagnostics.json").read_text(encoding="utf-8"))
    assert diagnostics["strict_clean"] is True
    assert diagnostics["counts"] == diagnostics["targets"]
    assert "cartesian_mixed_projection" not in diagnostics["sources"]
    assert not any(source.startswith("global_") for source in diagnostics["sources"])
