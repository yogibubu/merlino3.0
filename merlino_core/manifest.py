from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


def sha256_file(path: Path) -> str:
    """Return the SHA256 checksum of a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(path: Path, data: Mapping[str, Any]) -> Path:
    """Write a deterministic JSON manifest and return its path."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(dict(data), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return target
