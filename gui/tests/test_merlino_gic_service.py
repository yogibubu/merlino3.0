from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np
import pytest

from merlino_core.cli import build_parser as merlino_parser
from merlino_gic import (
    GICDefinition,
    GICDefinitionError,
    compare_gic_b_matrix_to_fortran,
    define_gics_from_cartesian,
    evaluate_gic_definition,
    read_gicforge_b_matrix,
    run_gicforge_python_fortran_contract,
    run_gicforge,
)
from merlino_gic.model import parse_gicforge_line
from merlino_fortran import resolve_backend
from merlino_gic.gic_symmetry import write_gic_symmetry_files
from merlino_fit.survibfit.primitives import Primitive
from merlino_fit.survibfit.cli import _python_local_gic_allowed


def test_gic_contract_cli_parser_accepts_contract_arguments():
    args = merlino_parser().parse_args(
        [
            "gic-contract",
            "--geometry",
            "water.xyz",
            "--workdir",
            "contract",
            "--json-out",
            "contract.json",
        ]
    )

    assert args.command == "gic-contract"
    assert str(args.geometry) == "water.xyz"
    assert str(args.workdir) == "contract"
    assert str(args.json_out) == "contract.json"


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
    assert manifest["backend"]["symmetrize"] is True


def test_run_gicforge_can_skip_symmetry_postcheck_and_removes_stale_outputs(tmp_path):
    executable = tmp_path / "fake_gicforge.sh"
    executable.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -e",
                "printf 'legacy report\\n' > provout",
                "printf 'raw gaussian input\\n' > gauin",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    (tmp_path / "provin").write_text("input\n", encoding="utf-8")
    (tmp_path / "gauin.symm").write_text("stale\n", encoding="utf-8")
    (tmp_path / "gicsym").write_text("stale\n", encoding="utf-8")

    result = run_gicforge(tmp_path, executable=executable, symmetrize=False)
    manifest = json.loads(result.manifest.read_text(encoding="utf-8"))

    assert "gauin" in result.files
    assert "gauin.symm" not in result.files
    assert "gicsym" not in result.files
    assert not (tmp_path / "gauin.symm").exists()
    assert manifest["backend"]["symmetrize"] is False


def test_gic_definition_can_be_reused_to_build_b_matrix_on_new_geometry(tmp_path):
    executable = tmp_path / "fake_gicforge.sh"
    executable.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -e",
                "printf ' Point Group from symm.f: C2v\\n' > provout",
                "cat > gauin.symm <<'EOF'",
                " A1Str0001=[ 0.70710678*R(  1,  2)+0.70710678*R(  1,  3)]",
                " B2Str0001=[ 0.70710678*R(  1,  2)-0.70710678*R(  1,  3)]",
                " A1Ang0001=[ 1.00000000*A(  2,  1,  3)]",
                "EOF",
                "cat > gicsym <<'EOF'",
                "name,irrep",
                "A1Str0001,A1",
                "B2Str0001,B2",
                "A1Ang0001,A1",
                "EOF",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    atoms = ("O", "H", "H")
    coords = np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [0.7570, 0.0000, 0.5860],
            [-0.7570, 0.0000, 0.5860],
        ],
        dtype=float,
    )

    definition = define_gics_from_cartesian(atoms, coords, workdir=tmp_path / "define", executable=executable)
    schema = definition.write(tmp_path / "gic_definition.json")
    restored = GICDefinition.read(schema)
    moved = coords.copy()
    moved[1, 0] += 0.01
    evaluation = evaluate_gic_definition(restored, moved, atomic_numbers=(8, 2, 2))

    assert restored.point_group == "C2v"
    assert restored.symmetrized is True
    assert restored.u_matrix.shape == (3, 3)
    assert len(restored.labels) == 3
    assert evaluation.values.shape == (3,)
    assert evaluation.b_matrix.shape == (3, 9)
    assert evaluation.irreps == ("A1", "B2", "A1")
    assert evaluation.point_group == "C2v"
    assert evaluation.symmetrized is True
    assert np.isfinite(evaluation.b_matrix).all()


def test_gic_definition_uses_backend_gauin_as_single_source(tmp_path):
    executable = tmp_path / "fake_gicforge.sh"
    executable.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -e",
                "printf ' Point Group from symm.f: C1\\n' > provout",
                "cat > gauin <<'EOF'",
                "0 1",
                "",
                " BackendOnly = R(  1,  3)",
                "EOF",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    atoms = ("O", "H", "H")
    coords = np.array(
        [
            [0.0000, 0.0000, 0.0000],
            [0.7570, 0.0000, 0.5860],
            [-0.7570, 0.0000, 0.5860],
        ],
        dtype=float,
    )

    definition = define_gics_from_cartesian(
        atoms,
        coords,
        workdir=tmp_path / "define",
        executable=executable,
        symmetrize=False,
    )
    evaluation = evaluate_gic_definition(definition, coords)

    assert definition.source == "gicforge"
    assert definition.symmetrized is False
    assert definition.point_group == "C1"
    assert definition.names == ("BackendOnly",)
    assert definition.provenance["backend"] == "gicforge"
    assert "xyzin_sha256" in definition.provenance
    assert "provin_sha256" in definition.provenance
    assert "gauin_sha256" in definition.provenance
    assert "backend_executable_sha256" in definition.provenance
    assert len(definition.primitives) == 1
    assert definition.primitives[0].kind == "bond"
    assert definition.primitives[0].atoms == (0, 2)
    assert definition.u_matrix.shape == (1, 1)
    assert not (tmp_path / "define" / "gauin.symm").exists()
    assert np.isclose(evaluation.values[0], np.linalg.norm(coords[0] - coords[2]))


def test_gic_definition_schema_validation_rejects_malformed_columns():
    bad = GICDefinition(
        atom_symbols=("H", "H"),
        atomic_numbers=(1, 1),
        reference_coordinates_angstrom=((0.0, 0.0, 0.0), (0.7, 0.0, 0.0)),
        primitives=(Primitive("bond", (0, 1)),),
        u_matrix=np.zeros((1, 1)),
        labels=("GIC001 zero",),
        names=("Zero",),
        irreps=("UNK",),
        symmetrized=False,
    )

    with pytest.raises(GICDefinitionError, match="zero-norm"):
        bad.write(Path("/tmp/should_not_write_gic_definition.json"))


def test_gicforge_b_matrix_triplets_match_python_evaluation(tmp_path):
    coords = np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], dtype=float)
    definition = GICDefinition(
        atom_symbols=("H", "H"),
        atomic_numbers=(1, 1),
        reference_coordinates_angstrom=tuple(tuple(row) for row in coords),
        primitives=(Primitive("bond", (0, 1)),),
        u_matrix=np.eye(1),
        labels=("GIC001 R(1,2)",),
        names=("BackendBond",),
        irreps=("UNK",),
        symmetrized=False,
    )
    python_b = evaluate_gic_definition(definition, coords).b_matrix
    bmat = tmp_path / "bmat.out"
    rows = ["# merlino.gicforge.bmatrix.v1", f"{python_b.shape[0]} {python_b.shape[1]}"]
    for row in range(python_b.shape[0]):
        for col in range(python_b.shape[1]):
            rows.append(f"{row + 1:8d}{col + 1:8d} {python_b[row, col]: .16E}".replace("E", "D"))
    bmat.write_text("\n".join(rows) + "\n", encoding="utf-8")

    parsed = read_gicforge_b_matrix(bmat)
    comparison = compare_gic_b_matrix_to_fortran(definition, coords, bmat)

    assert np.allclose(parsed, python_b)
    assert comparison.passed is True
    assert comparison.max_abs_diff == 0.0
    assert comparison.python_shape == python_b.shape
    assert comparison.fortran_shape == python_b.shape


def test_gicforge_python_fortran_contract_runs_real_backend(tmp_path):
    try:
        executable = resolve_backend("gicforge")
    except Exception as exc:
        pytest.skip(f"GICForge backend not available: {exc}")
    atoms = ("O", "H", "H")
    coords = np.array(
        [
            [0.000000, 0.000000, 0.000000],
            [0.757000, 0.000000, 0.586000],
            [-0.757000, 0.000000, 0.586000],
        ],
        dtype=float,
    )

    contract = run_gicforge_python_fortran_contract(
        atoms,
        coords,
        workdir=tmp_path / "contract",
        executable=executable,
    )

    assert contract.passed is True
    assert contract.contract_errors == ()
    assert contract.raw_b_matrix.passed is True
    assert contract.raw_b_matrix.max_abs_diff < 1.0e-7
    assert contract.raw_point_group == "C2v"
    assert contract.point_group == "C2v"
    assert "A1" in contract.irreps
    assert contract.raw_gic_count == 3
    assert contract.sym_gic_count == 3
    assert contract.totally_symmetric_count == 2
    assert len(contract.raw_names) == contract.raw_gic_count
    assert len(contract.sym_names) == contract.sym_gic_count
    assert len(contract.raw_labels) == contract.raw_gic_count
    assert len(contract.sym_labels) == contract.sym_gic_count
    assert all(signature.startswith(("bond:", "angle:")) for signature in contract.raw_primitive_signatures)
    assert sorted(contract.raw_primitive_signatures) == sorted(contract.sym_primitive_signatures)
    assert contract.raw_coordinate_kind_counts == {"angle": 1, "bond": 2}
    assert contract.sym_coordinate_kind_counts == contract.raw_coordinate_kind_counts


def test_gicforge_provout_final_summary_includes_ring_dihedral_coordinates(tmp_path):
    try:
        executable = resolve_backend("gicforge")
    except Exception as exc:
        pytest.skip(f"GICForge backend not available: {exc}")
    repo = Path(__file__).resolve().parents[2]
    xyz = repo / "benchmarks/semiexp_msr/inputs/cyclopentadiene/parent.xyz"
    lines = xyz.read_text(encoding="utf-8").splitlines()
    natoms = int(lines[0].strip())
    atoms: list[str] = []
    coords: list[list[float]] = []
    for raw in lines[2 : 2 + natoms]:
        fields = raw.split()
        atoms.append(fields[0])
        coords.append([float(value) for value in fields[1:4]])

    contract = run_gicforge_python_fortran_contract(
        atoms,
        np.asarray(coords, dtype=float),
        workdir=tmp_path / "contract",
        executable=executable,
    )
    provout = (tmp_path / "contract" / "raw" / "provout").read_text(encoding="utf-8", errors="replace")
    final_summary = provout[provout.index("Final GIC summary") :]

    assert contract.passed is True
    assert "Endocyclic Valence Angles" in provout
    assert "Exocyclic Dihedral Angles:" in provout
    assert "Endocyclic Dihedral Angles" in provout
    assert "RDef000" in final_summary
    assert "RPck000" in final_summary
    assert "QPck000" in final_summary
    assert provout.index("Endocyclic Valence Angles") < provout.index("Endocyclic Dihedral Angles")
    assert final_summary.index("RPck000") < final_summary.index("QPck000")


def test_gicforge_provout_reports_exocyclic_dihedral_count_before_butterfly(tmp_path):
    try:
        executable = resolve_backend("gicforge")
    except Exception as exc:
        pytest.skip(f"GICForge backend not available: {exc}")
    repo = Path(__file__).resolve().parents[2]
    xyz = repo / "merlino_fit/tests/data/polycyclics/saccharine.xyz"
    lines = xyz.read_text(encoding="utf-8").splitlines()
    natoms = int(lines[0].strip())
    atoms: list[str] = []
    coords: list[list[float]] = []
    for raw in lines[2 : 2 + natoms]:
        fields = raw.split()
        atoms.append(fields[0])
        coords.append([float(value) for value in fields[1:4]])

    contract = run_gicforge_python_fortran_contract(
        atoms,
        np.asarray(coords, dtype=float),
        workdir=tmp_path / "contract",
        executable=executable,
    )
    provout = (tmp_path / "contract" / "raw" / "provout").read_text(encoding="utf-8", errors="replace")
    provout_lines = provout.splitlines()
    butterfly_line = next(index for index, line in enumerate(provout_lines) if "Butterfly GNIC Around Bond" in line)

    assert contract.passed is True
    assert "Exocyclic Dihedral Angles:    0" in provout
    assert provout.index("Exocyclic Dihedral Angles:") < provout.index("Endocyclic Dihedral Angles")
    assert provout.index("Endocyclic Dihedral Angles") < provout.index("Butterfly GNIC Around Bond")
    assert provout_lines[butterfly_line - 1].strip()


def test_gicforge_parser_distinguishes_improper_dihedral_and_out_of_plane():
    impd = parse_gicforge_line(" ImpD0001 = D(  1,  2,  6,  4)")
    oupl = parse_gicforge_line(" OuPl0001 = U(  1,  2,  6,  4)")

    assert impd is not None
    assert oupl is not None
    assert impd[1] == [(1.0, Primitive("dihedral", (0, 1, 5, 3)))]
    assert oupl[1] == [(1.0, Primitive("out_of_plane", (0, 1, 5, 3)))]


def test_python_local_gic_requires_explicit_environment(monkeypatch):
    monkeypatch.delenv("MERLINO_ALLOW_PYTHON_LOCAL_GIC", raising=False)
    assert _python_local_gic_allowed() is False

    monkeypatch.setenv("MERLINO_ALLOW_PYTHON_LOCAL_GIC", "1")
    assert _python_local_gic_allowed() is True


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


def test_gic_symmetry_is_invariant_to_rigid_motion_and_equivalent_atom_swap(tmp_path):
    base = np.array(
        [
            [0.000000, 0.000000, 0.000000],
            [0.757000, 0.000000, 0.586000],
            [-0.757000, 0.000000, 0.586000],
        ],
        dtype=float,
    )
    theta = 0.37
    rotation = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    variants = {
        "base": (base, ("O", "H", "H"), ("R(  1,  2)", "R(  1,  3)", "A(  2,  1,  3)")),
        "moved": (base @ rotation.T + np.array([2.0, -1.0, 0.4]), ("O", "H", "H"), ("R(  1,  2)", "R(  1,  3)", "A(  2,  1,  3)")),
        "swapped": (base[[0, 2, 1]], ("O", "H", "H"), ("R(  1,  3)", "R(  1,  2)", "A(  3,  1,  2)")),
    }
    signatures = {}
    for name, (coords, atoms, primitive_texts) in variants.items():
        workdir = tmp_path / name
        workdir.mkdir()
        (workdir / "xyzin").write_text(
            "\n".join(
                [
                    "3",
                    name,
                    *(f"{atom} {x:.8f} {y:.8f} {z:.8f}" for atom, (x, y, z) in zip(atoms, coords)),
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        r12, r13, angle = primitive_texts
        (workdir / "gauin").write_text(
            "\n".join(
                [
                    "%chk=water.chk",
                    "",
                    "0 1",
                    "",
                    f" S1=[ 0.70710678*{r12}+0.70710678*{r13}]",
                    f" S2=[ 0.70710678*{r12}-0.70710678*{r13}]",
                    f" A1=[ 1.00000000*{angle}]",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        write_gic_symmetry_files(workdir)
        diagnostics = json.loads((workdir / "gic_symmetry_diagnostics.json").read_text(encoding="utf-8"))
        gicsym_rows = [
            tuple(line.split(",")[:2])
            for line in (workdir / "gicsym").read_text(encoding="utf-8").splitlines()[1:]
            if line.strip()
        ]
        signatures[name] = {
            "counts": diagnostics["counts"],
            "targets": diagnostics["targets"],
            "class_counts": diagnostics["class_counts"],
            "rows": tuple((row[0][:2], row[1]) for row in gicsym_rows),
        }

    assert signatures["base"]["counts"] == signatures["moved"]["counts"] == signatures["swapped"]["counts"]
    assert signatures["base"]["targets"] == signatures["moved"]["targets"] == signatures["swapped"]["targets"]
    assert signatures["base"]["class_counts"] == signatures["moved"]["class_counts"] == signatures["swapped"]["class_counts"]
    assert signatures["base"]["rows"] == signatures["moved"]["rows"] == signatures["swapped"]["rows"]
