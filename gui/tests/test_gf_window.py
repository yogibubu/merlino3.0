from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from advanced.gf_window import GFWindow
from merlino_fit.survibfit.primitives import Primitive
from merlino_gf import BOHR_TO_ANGSTROM
from merlino_gic import GICDefinition
from merlino_vpt2_vci.gaussian_qff import hessian_input_from_gaussian_fchk


@pytest.mark.usefixtures("qtbot")
def test_gf_window_defaults(tmp_path, qtbot):
    window = GFWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)

    assert window.windowTitle() == "Merlino GF / PED"
    assert "gauin.fchk" in window.fchk_edit.text()
    assert window.gic_schema_edit.text() == ""
    assert window.gic_geometry_edit.text() == ""
    assert window.gic_scale_edit.text() == ""


@pytest.mark.usefixtures("qtbot")
def test_gf_window_selects_latest_fchk(tmp_path, qtbot):
    old_fchk = tmp_path / "old.fchk"
    new_fchk = tmp_path / "gauin.fchk"
    old_fchk.write_text("old\n", encoding="utf-8")
    new_fchk.write_text("new\n", encoding="utf-8")
    os.utime(old_fchk, (1, 1))
    os.utime(new_fchk, (2, 2))

    window = GFWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window._select_latest_fchk(show_message=False)

    assert window.fchk_edit.text() == str(new_fchk)


@pytest.mark.usefixtures("qtbot")
def test_gf_window_runs_generated_gic_gf_on_gaussian_fchk(tmp_path, qtbot):
    fchk = Path("gui/tests/gaussian/h2o.fchk").resolve()
    window = GFWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.fchk_edit.setText(str(fchk))

    window.run_gf(show_message=False)

    text = window.output_text.toPlainText()
    assert "GF/PED from Merlino non-redundant GICs" in text
    assert "2169.878" in text
    assert "GIC001" in text
    assert "PED (%)" in text
    manifest = json.loads((tmp_path / "gf_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "merlino.run.v1"
    assert manifest["workflow"] == "gf"

    csv_dir = tmp_path / "gf_csv"
    written = window.export_csvs(csv_dir, show_message=False)
    assert (csv_dir / "gf_frequencies.csv").exists()
    assert "ped.csv" in written


@pytest.mark.usefixtures("qtbot")
def test_gf_window_runs_frozen_gic_definition_and_scaling(tmp_path, qtbot):
    fchk = Path("gui/tests/gaussian/h2o.fchk").resolve()
    canonical = hessian_input_from_gaussian_fchk(fchk)
    definition = GICDefinition(
        atom_symbols=("H", "O", "H"),
        atomic_numbers=(1, 8, 1),
        reference_coordinates_angstrom=tuple(tuple(row) for row in canonical.cartesian_coordinates_bohr * BOHR_TO_ANGSTROM),
        primitives=(
            Primitive("bond", (0, 1)),
            Primitive("bond", (1, 2)),
            Primitive("angle", (0, 1, 2)),
        ),
        u_matrix=np.eye(3),
        labels=("GIC001 R(1,2)", "GIC002 R(2,3)", "GIC003 A(1,2,3)"),
        names=("R1", "R2", "A1"),
        irreps=("A1", "A1", "A1"),
        point_group="C2v",
    )
    schema = definition.write(tmp_path / "gic_definition.json")
    scale = tmp_path / "scale.txt"
    scale.write_text("GIC003 0.90\n", encoding="utf-8")
    window = GFWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.fchk_edit.setText(str(fchk))
    window.gic_schema_edit.setText(str(schema))
    window.gic_scale_edit.setText(str(scale))

    window.run_gf(show_message=False)

    text = window.output_text.toPlainText()
    assert "Frozen GIC definition" in text
    assert "Pulay Hessian scaling: applied" in text
    manifest = json.loads((tmp_path / "gic_gf_manifest.json").read_text(encoding="utf-8"))
    assert manifest["workflow"] == "gic_gf"
    assert manifest["backend"]["coordinate_model"] == "frozen-gic-definition"


@pytest.mark.usefixtures("qtbot")
def test_gf_window_exports_report_and_roundtrips_preset(tmp_path, qtbot):
    window = GFWindow(tmp_path, tmp_path)
    qtbot.addWidget(window)
    window.output_text.setPlainText("sample report")
    report_path = tmp_path / "report.txt"

    written = window.export_report(report_path, show_message=False)

    assert written == report_path
    assert report_path.read_text(encoding="utf-8") == "sample report\n"

    window.fchk_edit.setText("/tmp/demo.fchk")
    window.gic_schema_edit.setText("/tmp/gic_definition.json")
    preset_path = tmp_path / "preset.json"

    saved = window.save_preset(preset_path, show_message=False)
    window.fchk_edit.setText("")
    window.gic_schema_edit.setText("")
    loaded = window.load_preset(preset_path, show_message=False)

    assert saved == preset_path
    assert loaded == preset_path
    assert window.fchk_edit.text() == "/tmp/demo.fchk"
    assert window.gic_schema_edit.text() == "/tmp/gic_definition.json"
