from __future__ import annotations

import json
from pathlib import Path

from merlino_gic import run_gicforge


def test_run_gicforge_collects_outputs_and_manifest(tmp_path):
    executable = tmp_path / "fake_gicforge.sh"
    executable.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -e",
                "printf 'legacy report\\n' > provout",
                "printf 'gaussian input\\n' > gauin",
                "printf 'vpt2 input\\n' > VPT2in",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    (tmp_path / "provin").write_text("input\n", encoding="utf-8")

    result = run_gicforge(tmp_path, executable=executable)

    assert result.files["provout"] == tmp_path / "provout"
    assert result.files["gicforge.out"].read_text(encoding="utf-8") == "legacy report\n"
    assert result.files["gauin"].read_text(encoding="utf-8") == "gaussian input\n"
    manifest = json.loads(result.manifest.read_text(encoding="utf-8"))
    assert manifest["workflow"] == "gicforge"
    assert manifest["outputs"]["gauin"] == str(tmp_path / "gauin")
    assert "provin" in manifest["inputs"]
