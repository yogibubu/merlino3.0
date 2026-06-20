#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python}"
"$PYTHON" -m pytest \
  gui/tests/test_merlino4_core_fortran.py \
  -q

