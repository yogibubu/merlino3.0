#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

echo "[freeze-check] root: $ROOT_DIR"
echo "[freeze-check] python: $(python3 --version)"
echo "[freeze-check] pytest: $(pytest --version)"
echo "[freeze-check] date: $(date '+%Y-%m-%d %H:%M:%S %Z')"

echo "[freeze-check] running GUI tests..."
QT_QPA_PLATFORM=offscreen pytest -q gui/tests

echo "[freeze-check] running merlino_fit tests..."
PYTHONPATH="$ROOT_DIR/merlino_fit" pytest -q merlino_fit/tests

echo "[freeze-check] OK"
