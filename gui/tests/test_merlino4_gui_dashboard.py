from __future__ import annotations

import pytest
from PySide6.QtCore import Qt

from merlino_gui import DashboardWindow, default_workflows
from merlino_gui.app import build_parser
from merlino_gui.dashboard import _semiexp_expert_diagnostics, workflow_detail_text, workflow_state_lines
from merlino_semiexp import read_observations


def test_default_workflows_include_new_scientific_areas():
    workflows = {item.workflow_id: item for item in default_workflows()}
    assert "vpt2_vci" in workflows
    assert "semiexp_geometry" in workflows
    assert workflows["vpt2_vci"].service == "merlino_vpt2_vci"
    assert workflows["gic_gf"].service == "merlino_gic + merlino_gf"
    assert workflows["semiexp_geometry"].service == "merlino_semiexp"
    assert "fortran77" in workflows["semiexp_geometry"].backends
    assert workflows["semiexp_geometry"].status == "standard solver"
    assert workflows["gic"].default_backend == "fortran77"
    assert workflows["gic"].status == "two-utility service available"
    assert workflows["gic_gf"].status == "frozen-GIC GF service available"


def test_workflow_detail_text_lists_contract_fields():
    workflow = default_workflows()[0]
    text = workflow_detail_text(workflow)
    assert workflow.title in text
    assert "Category:" in text
    assert "Inputs:" in text
    assert "Outputs:" in text
    assert workflow.service in text
    assert "Backend:" in text


def test_workflow_state_reports_manifest_outputs(tmp_path):
    run = tmp_path / "semiexp"
    run.mkdir()
    report = run / "semiexp_report.html"
    report.write_text("<html></html>\n", encoding="utf-8")
    manifest = run / "semiexp_manifest.json"
    manifest.write_text(
        """
{
  "workflow": "semiexperimental_geometry",
  "status": "completed",
  "run_dir": "%s",
  "outputs": {"html_report": "%s"}
}
"""
        % (run, report),
        encoding="utf-8",
    )
    workflow = {item.workflow_id: item for item in default_workflows()}["semiexp_geometry"]

    lines = workflow_state_lines(workflow, tmp_path)

    assert any("latest manifest: completed" in line for line in lines)
    assert any("outputs present: html_report" in line for line in lines)


def test_semiexp_expert_diagnostics_reads_manifest(tmp_path):
    manifest = tmp_path / "semiexp_manifest.json"
    manifest.write_text(
        """
{
  "parameters": {
    "convergence_reason": "objective_tolerance",
    "rank": 10,
    "incremental_rank": 10,
    "gicforge_calls": 3,
    "b_projector_secant_updates": 2
  },
  "outputs": {
    "influence": "/tmp/semiexp_influence.csv",
    "high_correlations": "/tmp/semiexp_high_correlations.csv"
  }
}
""",
        encoding="utf-8",
    )

    text = _semiexp_expert_diagnostics(manifest)

    assert "Expert diagnostics" in text
    assert "incremental_rank: 10" in text
    assert "b_projector_secant_updates: 2" in text
    assert "semiexp_high_correlations.csv" in text


@pytest.mark.usefixtures("qtbot")
def test_dashboard_lists_workflows(tmp_path, qtbot):
    window = DashboardWindow(tmp_path)
    qtbot.addWidget(window)

    assert window.windowTitle() == "Merlino 4.0"
    assert _workflow_tree_count(window) == len(default_workflows())
    listed = _workflow_tree_ids(window)
    assert "dvr" in listed
    assert "vpt2_vci" in listed
    assert "semiexp_geometry" in listed
    assert window.menuBar().actions()
    assert {window.menuBar().actions()[i].text() for i in range(window.menuBar().actions().__len__())} >= {
        "Project",
        "Structure",
        "Coordinates",
        "Vibrations",
        "Dynamics",
    }
    window.select_workflow("gic")
    assert not window.gic_panel.isHidden()
    assert window.semiexp_panel.isHidden()
    assert "gic-define" in window.gic_command.toPlainText()
    assert "gic-bmatrix" in window.gic_command.toPlainText()
    assert "gic-gf" in window.gic_command.toPlainText()
    assert window.gic_define_symmetry.isChecked()
    assert not window.gic_define_button.isEnabled()
    assert not window.gic_bmatrix_button.isEnabled()
    gic_xyz = tmp_path / "gic.xyz"
    gic_xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.957200",
                "H 0.926600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    window.gic_define_geometry.setText(str(gic_xyz))
    window.gic_define_schema.setText(str(tmp_path / "gic_definition.json"))
    window.gic_bmatrix_schema.setText(str(tmp_path / "gic_definition.json"))
    window.gic_bmatrix_geometry.setText(str(gic_xyz))
    window.gic_bmatrix_out.setText(str(tmp_path / "bmat.csv"))
    window.gic_gf_schema.setText(str(tmp_path / "gic_definition.json"))
    window.gic_gf_fchk.setText(str(tmp_path / "gauin.fchk"))
    window.gic_gf_report.setText(str(tmp_path / "gic_gf.txt"))
    window.gic_gf_csv_dir.setText(str(tmp_path / "gic_gf_csv"))
    assert window.gic_define_button.isEnabled()
    assert window.gic_bmatrix_button.isEnabled()
    assert window.gic_gf_button.isEnabled()
    assert window.gic_define_args()[:5] == [
        "gic-define",
        "--geometry",
        str(gic_xyz),
        "--out",
        str(tmp_path / "gic_definition.json"),
    ]
    assert window.gic_bmatrix_args()[:7] == [
        "gic-bmatrix",
        "--schema",
        str(tmp_path / "gic_definition.json"),
        "--geometry",
        str(gic_xyz),
        "--out",
        str(tmp_path / "bmat.csv"),
    ]
    assert window.gic_gf_args()[:5] == [
        "gic-gf",
        "--schema",
        str(tmp_path / "gic_definition.json"),
        "--fchk",
        str(tmp_path / "gauin.fchk"),
    ]
    window.select_workflow("gic_gf")
    assert not window.gic_panel.isHidden()
    assert window.semiexp_panel.isHidden()
    window.select_workflow("semiexp_geometry")
    assert window.backend_selector.isEnabled()
    assert not window.semiexp_panel.isHidden()
    assert window.gic_panel.isHidden()
    assert window.backend_selector.findText("python") >= 0
    assert window.backend_selector.findText("fortran77") >= 0
    assert window.semiexp_preview_table.columnCount() == 5
    window.backend_selector.setCurrentText("fortran77")
    assert window.selected_backends["semiexp_geometry"] == "fortran77"
    assert "selected: fortran77" in window.detail_view.toPlainText()
    parent_xyz = tmp_path / "parent.xyz"
    parent_xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.957200",
                "H 0.926600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    window.semiexp_xyz.setText(str(parent_xyz))
    window.semiexp_observations.setText(str(tmp_path / "isotopologues.toml"))
    window.semiexp_outdir.setText(str(tmp_path / "semiexp"))
    (tmp_path / "xyzin").write_text(parent_xyz.read_text(encoding="utf-8") + "\n#BASIC\ncharge 0\nmultiplicity 1\n", encoding="utf-8")
    window.semiexp_fixed.setText("GIC001")
    window.semiexp_gic_constraints.setText("QFIX=[GIC001+2*GIC002] Value=0.0; DR(Frozen,Value=0.0)=R[1,3]-R[1,2]")
    window.semiexp_fix_hydrogens.setChecked(True)
    window.semiexp_qm.setText("GIC002:1.0:0.1:qm")
    window.semiexp_classes.setText("CH:shared:R(1,2)|R(1,3);XYH:fixed:A(")
    args = window.semiexp_command_args()
    assert "--backend" in args
    assert args[args.index("--xyzin") + 1] == str(tmp_path / "xyzin")
    assert "fortran77" in args
    assert "--fix-hydrogens" in args
    assert args[args.index("--fixed") + 1] == "GIC001;QFIX=[GIC001+2*GIC002] Value=0.0;DR(Frozen,Value=0.0)=R[1,3]-R[1,2]"
    assert args[args.index("--prune-condition") + 1] == "0"
    assert args.count("--parameter-class") == 2
    assert window.semiexp_run_button.isEnabled()
    window.semiexp_iso_table.item(0, 2).setText("1000.0")
    window.semiexp_iso_table.item(0, 3).setText("800.0")
    window.semiexp_iso_table.item(0, 4).setText("600.0")
    toml = window.save_semiexp_observations_toml()
    assert toml.exists()
    assert "A_MHz = 1000.0" in toml.read_text(encoding="utf-8")
    xyzin_text = (tmp_path / "xyzin").read_text(encoding="utf-8")
    assert "#ISOTOPOLOGUES" in xyzin_text
    assert read_observations(tmp_path / "xyzin")[0].constants.A_MHz == 1000.0
    job = window.save_semiexp_job_toml()
    assert job is not None
    job_text = job.read_text(encoding="utf-8")
    assert 'schema = "merlino.semiexp.job.v1"' in job_text
    assert '["O", 0, 0, 0]' in job_text
    assert "fix_hydrogen_parameters = true" in job_text
    assert "gic_constraints = [" in job_text
    assert "QFIX=[GIC001+2*GIC002] Value=0.0" in job_text
    assert "[[isotopologues]]" in job_text
    assert "[isotopologues.definition]" in job_text
    assert "[files]" not in job_text
    preset = window.save_semiexp_preset()
    assert preset.exists()
    window.semiexp_fixed.clear()
    window.semiexp_gic_constraints.clear()
    window.semiexp_fix_hydrogens.setChecked(False)
    window.load_semiexp_preset(preset)
    assert window.semiexp_fixed.text() == "GIC001"
    assert window.semiexp_gic_constraints.text() == "QFIX=[GIC001+2*GIC002] Value=0.0; DR(Frozen,Value=0.0)=R[1,3]-R[1,2]"
    assert window.semiexp_fix_hydrogens.isChecked()


@pytest.mark.usefixtures("qtbot")
def test_semiexp_dashboard_preview_validate_and_conditioning(tmp_path, qtbot):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.957200",
                "H 0.926600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    obs = tmp_path / "isotopologues.toml"
    obs.write_text(
        "\n".join(
            [
                "[[isotopologues]]",
                'label = "parent"',
                "[isotopologues.constants]",
                "A_MHz = 822180.189172425",
                "B_MHz = 437776.592728178",
                "C_MHz = 285669.514220614",
                "[isotopologues.vibrational_correction]",
                "delta_A_MHz = 0.0",
                "delta_B_MHz = 0.0",
                "delta_C_MHz = 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    window = DashboardWindow(tmp_path)
    qtbot.addWidget(window)
    window.select_workflow("semiexp_geometry")
    window.semiexp_xyz.setText(str(xyz))
    window.semiexp_observations.setText(str(obs))
    window.semiexp_outdir.setText(str(tmp_path / "run"))
    (tmp_path / "xyzin").write_text(xyz.read_text(encoding="utf-8") + "\n#BASIC\ncharge 0\nmultiplicity 1\n", encoding="utf-8")

    window.preview_semiexp_gics()
    window.validate_semiexp_input()
    window.preview_semiexp_conditioning()
    job = window.save_semiexp_job_toml()

    text = window.semiexp_command.toPlainText()
    assert window.semiexp_preview_table.rowCount() > 0
    assert "input validation: OK" in text
    assert "condition number" in text
    assert job is not None
    assert "#ISOTOPOLOGUES" in (tmp_path / "xyzin").read_text(encoding="utf-8")
    assert read_observations(tmp_path / "xyzin")[0].constants.A_MHz == pytest.approx(822180.189172425)


@pytest.mark.usefixtures("qtbot")
def test_semiexp_dashboard_can_use_inline_isotopologue_job(tmp_path, qtbot):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.957200",
                "H 0.926600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    window = DashboardWindow(tmp_path)
    qtbot.addWidget(window)
    window.select_workflow("semiexp_geometry")
    window.semiexp_xyz.setText(str(xyz))
    window.semiexp_observations.clear()
    window.semiexp_outdir.setText(str(tmp_path / "run"))
    window.semiexp_iso_table.item(0, 2).setText("822180.189172425")
    window.semiexp_iso_table.item(0, 3).setText("437776.592728178")
    window.semiexp_iso_table.item(0, 4).setText("285669.514220614")

    args = window.semiexp_command_args()
    job = window.save_semiexp_job_toml()

    assert args[:3] == ["semiexp", "--job", str(tmp_path / "semiexp_job.mse.toml")]
    assert args[args.index("--xyzin") + 1] == str(tmp_path / "xyzin")
    assert "--observations" not in args
    assert job is not None
    text = job.read_text(encoding="utf-8")
    assert "[[isotopologues]]" in text
    assert "[isotopologues.definition]" in text
    assert "[files]" not in text


@pytest.mark.usefixtures("qtbot")
def test_semiexp_dashboard_builds_advanced_sefit_command_and_job(tmp_path, qtbot):
    xyz = tmp_path / "parent.xyz"
    xyz.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.000000 0.957200",
                "H 0.926600 0.000000 -0.239600",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    obs = tmp_path / "isotopologues.toml"
    obs.write_text(
        "\n".join(
            [
                "[[isotopologues]]",
                'label = "parent"',
                "[isotopologues.constants]",
                "A_MHz = 822180.189172425",
                "B_MHz = 437776.592728178",
                "C_MHz = 285669.514220614",
                "[isotopologues.vibrational_correction]",
                "delta_A_MHz = 0.0",
                "delta_B_MHz = 0.0",
                "delta_C_MHz = 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    window = DashboardWindow(tmp_path)
    qtbot.addWidget(window)
    window.select_workflow("semiexp_geometry")
    window.backend_selector.setCurrentText("fortran77")
    window.semiexp_xyz.setText(str(xyz))
    window.semiexp_observations.setText(str(obs))
    window.semiexp_outdir.setText(str(tmp_path / "run"))
    window.semiexp_coordinate_model.setCurrentText("cartesian_symmetry")
    window.semiexp_observable.setCurrentText("moments")
    window.semiexp_components.setCurrentText("AB")
    window.semiexp_prune_condition.setValue(1000.0)
    window.semiexp_robust_loss.setCurrentText("huber")
    window.semiexp_robust_scale.setValue(2.5)
    window.semiexp_leave_one_out.setChecked(True)
    window.semiexp_checkpoint.setText(str(tmp_path / "checkpoint.json"))
    window.semiexp_restart.setText(str(tmp_path / "restart.json"))
    window.semiexp_fix_hydrogens.setChecked(True)
    window.semiexp_gic_constraints.setText("DR(Frozen,Value=0.0)=R[1,3]-R[1,2]")
    window.semiexp_qm.setText("R(1,2):0.9572:0.002:BDPCS3")
    window.semiexp_classes.setText("OH:shared:R(1,2)|R(1,3)")

    args = window.semiexp_command_args()

    assert args[args.index("--backend") + 1] == "fortran77"
    assert args[args.index("--xyzin") + 1] == str(tmp_path / "xyzin")
    assert args[args.index("--coordinate-model") + 1] == "cartesian_symmetry"
    assert args[args.index("--observable") + 1] == "moments"
    assert args[args.index("--rotational-components") + 1] == "AB"
    assert args[args.index("--prune-condition") + 1] == "1000"
    assert args[args.index("--robust-loss") + 1] == "huber"
    assert args[args.index("--robust-scale") + 1] == "2.5"
    assert "--leave-one-out" in args
    assert args[args.index("--checkpoint") + 1] == str(tmp_path / "checkpoint.json")
    assert args[args.index("--restart") + 1] == str(tmp_path / "restart.json")
    assert "--fix-hydrogens" in args
    assert args[args.index("--fixed") + 1] == "DR(Frozen,Value=0.0)=R[1,3]-R[1,2]"
    assert args[args.index("--qm-predicate") + 1] == "R(1,2):0.9572:0.002:BDPCS3"
    assert args[args.index("--parameter-class") + 1] == "OH:shared:R(1,2)|R(1,3)"
    job = window.save_semiexp_job_toml()
    assert job is not None
    text = job.read_text(encoding="utf-8")
    assert 'backend = "fortran77"' in text
    assert 'coordinate_model = "cartesian_symmetry"' in text
    assert 'rotational_components = "AB"' in text
    assert 'robust_loss = "huber"' in text
    assert "robust_scale = 2.5" in text
    assert "leave_one_out = true" in text
    assert f'checkpoint = "{tmp_path / "checkpoint.json"}"' in text
    assert f'restart = "{tmp_path / "restart.json"}"' in text
    assert "fix_hydrogen_parameters = true" in text
    assert "gic_constraints = [" in text
    assert 'pattern = "R(1,2)"' in text


def test_dashboard_launcher_parser_accepts_workdir(tmp_path):
    args = build_parser().parse_args(["--workdir", str(tmp_path)])
    assert args.workdir == tmp_path


def _workflow_tree_count(window):
    total = 0
    for top_idx in range(window.workflow_list.topLevelItemCount()):
        total += window.workflow_list.topLevelItem(top_idx).childCount()
    return total


def _workflow_tree_ids(window):
    ids = set()
    for top_idx in range(window.workflow_list.topLevelItemCount()):
        parent = window.workflow_list.topLevelItem(top_idx)
        for child_idx in range(parent.childCount()):
            ids.add(parent.child(child_idx).data(0, Qt.UserRole))
    return ids
