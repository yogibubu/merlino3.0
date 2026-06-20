#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python}"
"$PYTHON" -m pytest \
  gui/tests/test_merlino4_gui_dashboard.py \
  gui/tests/test_vpt2_vci_window.py \
  gui/tests/test_dvr_window.py \
  -q
