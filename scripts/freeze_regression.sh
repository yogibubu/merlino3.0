#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python}"
"$PYTHON" -m pytest \
  gui/tests/test_scientific_contracts.py \
  gui/tests/test_merlino4_core_infrastructure.py \
  gui/tests/test_merlino4_gui_dashboard.py \
  -q

