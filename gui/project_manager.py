from pathlib import Path

from geometry.rotational_pipeline import rotational_pipeline
from geometry.thermo_pipeline import thermo_pipeline
from geometry.rovib_pipeline import rovib_pipeline
from merlino_fit.topology.test_topology import run_topology_on_xyzin


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

    def run_topology(self) -> None:
        self.validate_xyzin()
        self.topology_result = run_topology_on_xyzin(str(self.xyzin))

    def run_rovib(self) -> None:
        self.validate_xyzin()
        opts = dict(self.rovib_options or {})
        self.rovib_result = rovib_pipeline(str(self.xyzin), **opts)

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
