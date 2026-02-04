from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import json


@dataclass
class GuiSettings:
    dos_do_vib: bool = True
    dos_do_rovib: bool = True
    dos_emin: float = 0.0
    dos_emax: float = 8000.0
    dos_bin: float = 50.0
    dos_vmax: int = 6
    dos_ncap: str = "10"
    dos_T: float = 298.15
    dos_emax_rot: float | None = None
    dos_jmax: int | None = None

    @classmethod
    def from_json(cls, path: Path, defaults: "GuiSettings" | None = None) -> "GuiSettings":
        if defaults is None:
            defaults = cls()
        if not path.exists():
            return defaults
        try:
            data = json.loads(path.read_text())
        except Exception:
            return defaults
        values = asdict(defaults)
        for key, val in data.items():
            if key in values:
                values[key] = val
        return cls(**values)

    def to_json(self, path: Path):
        path.write_text(json.dumps(asdict(self), indent=2))
