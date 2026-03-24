from pathlib import Path

from geometry.rotational_pipeline import rotational_pipeline
from geometry.vibrational import vib_from_xyzin
from geometry.thermo_vib import _read_blocks, _parse_vibrational_block
from geometry.rotational_pipeline import _write_vibrational_report
from geometry.coriolis import run_coriolis_from_vibin
from geometry.thermo_pipeline import thermo_pipeline
from geometry.rovib_pipeline import rovib_pipeline
from topology.test_topology import run_topology_on_xyzin


class ProjectManager:
    def __init__(self, working):
        self.working = Path(working)
        self.xyzin = self.working / "xyzin"

        # store last results (optional, non-breaking)
        self.rotational_result = None
        self.thermo_result = None
        self.topology_result = None
        self.rovib_result = None
        self.rovib_options = {}

    # ==========================================================
    # Validation
    # ==========================================================
    def validate_xyzin(self) -> None:
        if not self.xyzin.exists():
            raise FileNotFoundError(f"xyzin not found: {self.xyzin}")

    # ==========================================================
    # Explicit pipelines (atomic)
    # ==========================================================
    def run_rotational(self) -> None:
        self.validate_xyzin()
        self.rotational_result = rotational_pipeline(str(self.xyzin))

    def run_thermo(self) -> None:
        self.validate_xyzin()
        self.thermo_result = thermo_pipeline(str(self.xyzin))

    def run_topology(
        self,
        *,
        symm_tol: float = 1.0e-3,
        symmetrize_coords: bool = False,
    ) -> None:
        self.validate_xyzin()
        self.topology_result = run_topology_on_xyzin(
            str(self.xyzin),
            symm_tol=float(symm_tol),
            symmetrize_coords=bool(symmetrize_coords),
        )

    def run_rovib(self) -> None:
        self.validate_xyzin()
        opts = dict(self.rovib_options or {})
        self.rovib_result = rovib_pipeline(str(self.xyzin), **opts)

    def run_vibrational(self) -> None:
        self.validate_xyzin()
        fchkin_path = self.xyzin.parent / "fchkin"
        if fchkin_path.exists():
            vib = vib_from_xyzin(str(self.xyzin))
            if vib is None:
                return
            vibin_path = self.xyzin.parent / "vibin"
            coriolis_entries = None
            if vibin_path.exists():
                try:
                    coriolis_entries = run_coriolis_from_vibin(
                        vibin_path=str(vibin_path),
                        k_threshold=0.1,
                        keep_hf_modes=True,
                    )
                except Exception:
                    coriolis_entries = None
            vib_report_path = self.xyzin.parent / "vibrational.report"
            _write_vibrational_report(vib_report_path, vib, coriolis_entries=coriolis_entries)
            self.rovib_result = vib
            return

        # Fallback: use #VIBRATIONAL block if present (e.g., from Gaussian log)
        try:
            _basic, vib_lines = _read_blocks(str(self.xyzin))
            if not vib_lines:
                return
            freq = _parse_vibrational_block(vib_lines)
        except Exception:
            return

        vib_report_path = self.xyzin.parent / "vibrational.report"
        lines = [
            "VIBRATIONAL REPORT",
            "------------------",
            "Source: #VIBRATIONAL block (no fchkin)",
            f"Modes: {len(freq)}",
            "",
            "Frequencies (cm^-1):",
        ]
        for i, f in enumerate(freq, start=1):
            lines.append(f"{i:4d}  {float(f): .6f}")
        vib_report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ==========================================================
    # Canonical workflow (ORDERED, SAFE)
    # ==========================================================
    def run_full_workflow(self) -> dict:
        """
        Canonical Merlino workflow:

        1) Rotational
        2) Thermo
        3) Topology
        """
        self.run_rotational()
        self.run_thermo()
        self.run_topology()

        return {
            "rotational": self.rotational_result,
            "thermo": self.thermo_result,
            "topology": self.topology_result,
        }

    # ==========================================================
    # Flexible execution (advanced / expert use)
    # ==========================================================
    def run(self, pipelines):
        """
        Run selected pipelines.

        NOTE:
        - Order is enforced only in run_full_workflow().
        - This method is for expert / partial workflows.
        """
        self.validate_xyzin()

        if isinstance(pipelines, str):
            pipelines = [pipelines]

        for p in pipelines:
            p = p.lower()
            if p == "rotational":
                self.run_rotational()
            elif p == "thermo":
                self.run_thermo()
            elif p == "topology":
                self.run_topology()
            elif p == "rovib":
                self.run_rovib()
            else:
                raise ValueError(f"Unknown pipeline: {p}")
