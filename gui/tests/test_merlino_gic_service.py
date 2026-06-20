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


def test_gic_symmetry_keeps_cs_mirror_with_identity_permutation(tmp_path):
    (tmp_path / "xyzin").write_text(
        "\n".join(
            [
                "3",
                "planar Cs asymmetric triatomic",
                "C 0.000000 0.000000 0.000000",
                "O 1.200000 0.000000 0.000000",
                "N 0.250000 1.100000 0.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "gauin").write_text(
        "\n".join(
            [
                "%chk=cs.chk",
                "",
                "0 1",
                "",
                " R1=[ 1.00000000*R(  1,  2)]",
                " R2=[ 1.00000000*R(  1,  3)]",
                " A1=[ 1.00000000*A(  2,  1,  3)]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    write_gic_symmetry_files(tmp_path)
    diagnostics = json.loads((tmp_path / "gic_symmetry_diagnostics.json").read_text(encoding="utf-8"))
    gicsym = (tmp_path / "gicsym").read_text(encoding="utf-8")

    assert diagnostics["operation_order"] == ["E", "sigma_xy"]
    assert diagnostics["targets"] == {"A'": 3, "A''": 0}
    assert diagnostics["counts"] == {"A'": 3, "A''": 0}
    assert "A'Str0001,A'" in gicsym
    assert "cartesian_mixed_projection" not in diagnostics["sources"]
    assert not any(source.startswith("global_") for source in diagnostics["sources"])


def test_gic_symmetry_preserves_cyclopentadiene_irrep_and_class_counts(tmp_path):
    (tmp_path / "xyzin").write_text(
        "\n".join(
            [
                "11",
                "cyclopentadiene MSR geometry",
                "C 0.0000000648 1.21099405 -0.0000437769",
                "C 1.16966652 0.34210493 0.0000360217",
                "C -1.16966649 0.34210505 0.0000319263",
                "C 0.731740874 -0.975616979 -0.0000177892",
                "C -0.731740976 -0.975616904 -0.0000203515",
                "H 2.26478781 0.374863776 0.0001274146",
                "H -2.26478777 0.374864010 0.0001200854",
                "H 1.26945754 -1.92001108 0.0000062677",
                "H -1.26945774 -1.92001095 0.0000018591",
                "H 0.0000585394 1.87871494 -0.8689120000",
                "H -0.0000583448 1.87871720 0.8688227060",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "gauin").write_text(
        "\n".join(
            [
                "%chk=gicforge.chk",
                "",
                "0 1",
                "",
                " Stre0001 = R(  1,  2)",
                " Stre0002 = R(  1,  3)",
                " Stre0003 = R(  2,  4)",
                " Stre0004 = R(  3,  5)",
                " Stre0005 = R(  4,  5)",
                " Stre0006 = R(  1, 10)",
                " Stre0007 = R(  1, 11)",
                " Stre0008 = R(  2,  6)",
                " Stre0009 = R(  3,  7)",
                " Stre0010 = R(  4,  8)",
                " Stre0011 = R(  5,  9)",
                " SymD0001 =[ 0.81650*A( 10,  1, 11)- 0.40825*A(  2,  1, 10)- 0.40825*A(  2,  1, 11)]",
                " Rock0002 =[ 0.70711*A(  2,  1, 10)- 0.70711*A(  2,  1, 11)]",
                " SymD0003 =[ 0.81650*A( 10,  1, 11)- 0.40825*A(  3,  1, 10)- 0.40825*A(  3,  1, 11)]",
                " Rock0004 =[ 0.70711*A(  3,  1, 10)- 0.70711*A(  3,  1, 11)]",
                " Rock0005 =[ 0.70711*A(  6,  2,  4)- 0.70711*A(  6,  2,  1)]",
                " Rock0006 =[ 0.70711*A(  7,  3,  5)- 0.70711*A(  7,  3,  1)]",
                " Rock0007 =[ 0.70711*A(  8,  4,  5)- 0.70711*A(  8,  4,  2)]",
                " Rock0008 =[ 0.70711*A(  9,  5,  4)- 0.70711*A(  9,  5,  3)]",
                " RDef0009 =[ 0.63246*A(  5,  3,  1)- 0.51167*A(  3,  1,  2)+ 0.19544*A(  1,  2,  4)+ 0.19544*A(  2,  4,  5)- 0.51167*A(  4,  5,  3)]",
                " RDef0010 =[ 0.00000*A(  5,  3,  1)+ 0.37175*A(  3,  1,  2)- 0.60150*A(  1,  2,  4)+ 0.60150*A(  2,  4,  5)- 0.37175*A(  4,  5,  3)]",
                " RPck0001(Inactive)=[   0.63246*D(  1,  2,  4,  5)-0.51167*D(  2,  4,  5,  3)+0.19544*D(  4,  5,  3,  1)+0.19544*D(  5,  3,  1,  2)-0.51167*D(  3,  1,  2,  4)]",
                " RPck0002(Inactive)=[   0.00000*D(  1,  2,  4,  5)+0.37175*D(  2,  4,  5,  3)-0.60150*D(  4,  5,  3,  1)+0.60150*D(  5,  3,  1,  2)-0.37175*D(  3,  1,  2,  4)]",
                " ImpD0001 = D(  1,  2,  6,  4)",
                " ImpD0002 = D(  1,  3,  7,  5)",
                " ImpD0003 = D(  2,  4,  8,  5)",
                " ImpD0004 = D(  3,  5,  9,  4)",
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
    diagnostics = json.loads((tmp_path / "gic_symmetry_diagnostics.json").read_text(encoding="utf-8"))
    gicsym = (tmp_path / "gicsym").read_text(encoding="utf-8")

    assert first == second
    assert diagnostics["counts"] == {"A1": 10, "A2": 4, "B1": 5, "B2": 8}
    assert diagnostics["b_ranks"] == diagnostics["targets"]
    assert diagnostics["class_counts"] == diagnostics["class_targets"]
    assert diagnostics["class_counts"] == {"bond": 11, "angle": 10, "dihedral": 6}
    a1_lines = [line for line in gicsym.splitlines() if ",A1," in line]
    assert sum(line.startswith("A1Str") for line in a1_lines) > 0
    assert sum(line.startswith("A1Ang") for line in a1_lines) > 0
    assert "Tor" in gicsym
    assert "cartesian_mixed_projection" not in diagnostics["sources"]
    assert not any(source.startswith("global_") for source in diagnostics["sources"])
