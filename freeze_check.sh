#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-python}"
cd "$ROOT_DIR"

echo "[freeze-check] root: $ROOT_DIR"
echo "[freeze-check] python: $("$PYTHON" --version)"
echo "[freeze-check] pytest: $("$PYTHON" -m pytest --version)"
echo "[freeze-check] date: $(date '+%Y-%m-%d %H:%M:%S %Z')"

echo "[freeze-check] running GUI tests..."
QT_QPA_PLATFORM=offscreen "$PYTHON" -m pytest -q gui/tests

echo "[freeze-check] running merlino_fit tests..."
PYTHONPATH="$ROOT_DIR/merlino_fit" "$PYTHON" -m pytest -q merlino_fit/tests

echo "[freeze-check] OK"
