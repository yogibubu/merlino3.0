#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_NAME="$(basename "$ROOT_DIR")"
STAMP="$(date '+%Y%m%d_%H%M%S')"
PARENT_DIR="$(dirname "$ROOT_DIR")"
OUT_DIR="$ROOT_DIR/archives/freezes"

RUN_CHECKS=1
for arg in "$@"; do
  case "$arg" in
    --skip-check) RUN_CHECKS=0 ;;
    --out-dir=*) OUT_DIR="${arg#*=}" ;;
    *)
      echo "Unknown option: $arg"
      echo "Usage: $0 [--skip-check] [--out-dir=/path]"
      exit 2
      ;;
  esac
done

ARCHIVE_BASE="${PROJECT_NAME}_freeze_${STAMP}"
ARCHIVE_PATH="${OUT_DIR%/}/${ARCHIVE_BASE}.tar.gz"
SHA_PATH="${OUT_DIR%/}/${ARCHIVE_BASE}.sha256"

cd "$ROOT_DIR"

if [[ "$RUN_CHECKS" -eq 1 ]]; then
  echo "[freeze-pack] running freeze_check.sh..."
  ./freeze_check.sh
fi

echo "[freeze-pack] creating archive: $ARCHIVE_PATH"
tar \
  --exclude=".git" \
  --exclude=".pytest_cache" \
  --exclude="__pycache__" \
  --exclude="*/__pycache__" \
  --exclude=".DS_Store" \
  --exclude="work" \
  --exclude="working" \
  --exclude="archives/freezes" \
  --exclude="*_freeze_*.tar.gz" \
  --exclude="*_freeze_*.sha256" \
  -czf "$ARCHIVE_PATH" \
  -C "$PARENT_DIR" \
  "$PROJECT_NAME"

echo "[freeze-pack] writing checksum: $SHA_PATH"
shasum -a 256 "$ARCHIVE_PATH" > "$SHA_PATH"

echo "[freeze-pack] done"
echo "  archive:  $ARCHIVE_PATH"
echo "  checksum: $SHA_PATH"
