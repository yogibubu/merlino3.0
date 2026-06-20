#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python}"
"$PYTHON" -m pytest \
  gui/tests/test_scientific_contracts.py::test_semiexperimental_geometry_fit_reduces_rotational_residuals \
  gui/tests/test_scientific_contracts.py::test_semiexperimental_parameter_classes_share_and_fix_parameters \
  gui/tests/test_scientific_contracts.py::test_semiexperimental_gic_preview_and_html_report \
  gui/tests/test_merlino4_core_infrastructure.py::test_merlino_cli_semiexp \
  -q

