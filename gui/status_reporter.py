from __future__ import annotations

from pathlib import Path
from PySide6.QtCore import QUrl


class StatusReporter:
    def __init__(self, working_dir: Path):
        self.working_dir = working_dir

    def update_status_label(
        self,
        label,
        input_type,
        source,
        vib_q,
        rovib_q,
        dos_emin,
        dos_emax,
        dos_bin,
        dos_T,
        error,
    ):
        parts = []
        parts.append(f"Input: {input_type or 'unknown'}")
        if source:
            parts.append(f"Source: {source}")
        if vib_q is not None:
            parts.append(f"Q_vib: {vib_q:.6e}")
        if rovib_q is not None:
            parts.append(f"Q_rovib: {rovib_q:.6e}")
        parts.append(f"Emin/Emax/bin: {dos_emin:.1f}/{dos_emax:.1f}/{dos_bin:.1f} cm-1")
        parts.append(f"T: {dos_T:.2f} K")
        files = []
        for name in (
            "dos_vib.dat",
            "vib_qt.dat",
            "dos_rovib.dat",
            "rovib_qt.dat",
            "n_vib_ts.dat",
            "n_rovib_ts.dat",
            "gui.log",
            "summary.txt",
        ):
            fpath = self.working_dir / name
            if fpath.exists():
                url = QUrl.fromLocalFile(str(self.working_dir)).toString()
                files.append(f"<a href=\"{url}\">{name}</a>")
        if files:
            parts.append("Files: " + " | ".join(files))
        if error:
            parts.append(f"Warning: {error}")
        line1 = " | ".join(parts[:4])
        line2 = " | ".join(parts[4:]) if len(parts) > 4 else ""
        if error:
            line1 = "⚠️ " + line1
        text = line1 if not line2 else (line1 + "<br>" + line2)
        label.setText(text)

    def write_summary(
        self,
        input_type,
        source,
        dos_emin,
        dos_emax,
        dos_bin,
        dos_vmax,
        dos_ncap,
        dos_T,
        vib_q,
        rovib_q,
        error,
    ):
        try:
            lines = []
            lines.append("MERLINO GUI SUMMARY\n")
            lines.append(f"Input type: {input_type or 'unknown'}\n")
            lines.append(f"Source: {source}\n")
            lines.append(f"Emin_cm1: {dos_emin}\n")
            lines.append(f"Emax_cm1: {dos_emax}\n")
            lines.append(f"Bin_cm1: {dos_bin}\n")
            lines.append(f"vmax: {dos_vmax}\n")
            lines.append(f"ncap: {dos_ncap}\n")
            lines.append(f"T_K: {dos_T}\n")
            if vib_q is not None:
                lines.append(f"Q_vib: {vib_q:.12e}\n")
            if rovib_q is not None:
                lines.append(f"Q_rovib: {rovib_q:.12e}\n")
            if error:
                lines.append(f"Warning: {error}\n")
            files = []
            for name in (
                "dos_vib.dat",
                "vib_qt.dat",
                "dos_rovib.dat",
                "rovib_qt.dat",
                "n_vib_ts.dat",
                "n_rovib_ts.dat",
                "gui.log",
            ):
                if (self.working_dir / name).exists():
                    files.append(name)
            if files:
                lines.append("Files:\n")
                for f in files:
                    lines.append(f"- {f}\n")
            (self.working_dir / "summary.txt").write_text("".join(lines), encoding="utf-8")
        except Exception:
            pass
