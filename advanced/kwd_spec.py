"""
gui_new.kwd_spec

Keyword specification for the Merlino 3.0 GUI.

This module defines the controlled vocabulary used by the GUI:
- available keywords
- semantic grouping
- mutual exclusivity
- combinability

NO execution logic.
NO validation.
NO workflow decisions.
"""

KWD_SPEC = {

    # =========================================================
    # INPUT
    # =========================================================
    "input": {
        "panel": "Input",
        "description": "Molecular input and basic molecular information",
        "keywords": {
            "INPUT_FORMAT": {
                "values": [
                    "xyz",
                    "xyz_enriched",
                    "gaussian_output",
                    "molpro_output",
                    "mrcc_output",
                    "smiles",
                ],
                "exclusive": True,
                "description": "Type of molecular input",
            },
            "CHARGE": {
                "type": "int",
                "default": 0,
                "description": "Total molecular charge",
            },
            "MULTIPLICITY": {
                "type": "int",
                "default": 1,
                "description": "Spin multiplicity",
            },
        },
    },

    # =========================================================
    # GEOMETRY / TOPOLOGY PREPROCESSING
    # =========================================================
    "geometry_topology": {
        "panel": "Geometry / Topology",
        "description": "Preprocessing of molecular structure and connectivity",
        "keywords": {
            "FINDFR": {
                "type": "flag",
                "default": False,
                "description": "Identify molecular fragments",
            },
            "JOINFR": {
                "type": "flag",
                "default": False,
                "description": "Join molecular fragments",
            },
            "HBOND": {
                "type": "flag",
                "default": False,
                "description": "Detect/report non-covalent H-bond targets without adding them to GIC topology",
            },
        },
    },

    # =========================================================
    # GNIC
    # =========================================================
    "gnic": {
        "panel": "GNIC",
        "description": "Geometrical preprocessing and internal coordinate handling",
        "keywords": {
            "GNIC_MODE": {
                "values": [
                    None,
                    "NoOneDih",
                    "BDPCS3",
                ],
                "exclusive": True,
                "default": None,
                "description": "GNIC preprocessing mode",
            },
        },
    },

    # =========================================================
    # ELECTRONIC STRUCTURE MODEL
    # =========================================================
    "method": {
        "panel": "Method",
        "description": "Electronic structure model",
        "keywords": {
            "METHOD": {
                "values": [
                    "PCS0",
                    "PCS1",
                    "HPCS2",
                    "DPCS3",
                ],
                "exclusive": True,
                "description": "Composite electronic structure scheme",
            },
        },
    },

    # =========================================================
    # CALCULATION TYPE
    # =========================================================
    "calculation": {
        "panel": "Calculation",
        "description": "Requested computational steps",
        "keywords": {
            "SP": {
                "type": "flag",
                "default": False,
                "description": "Single-point energy calculation",
            },
            "OPT": {
                "type": "flag",
                "default": False,
                "description": "Geometry optimization",
            },
            "HARM": {
                "type": "flag",
                "default": False,
                "description": "Harmonic vibrational analysis",
            },
            "ANH": {
                "type": "flag",
                "default": False,
                "description": "Anharmonic vibrational analysis",
            },
        },
    },

    # =========================================================
    # SYMMETRY
    # =========================================================
    "symmetry": {
        "panel": "Symmetry",
        "description": "Molecular symmetry analysis options",
        "keywords": {
            "GICSYM": {
                "type": "flag",
                "default": False,
                "description": "Symmetrize GIC blocks for downstream modules",
            },
            "SYCART": {
                "type": "flag",
                "default": False,
                "description": "Write symmetrized Cartesian coordinates without changing the input frame",
            },
            "SYMMALL": {
                "type": "flag",
                "default": False,
                "description": "Legacy alias for GIC block symmetrization",
            },
            "LOOSE": {
                "type": "flag",
                "default": False,
                "description": "Loose tolerances for symmetry detection",
            },
            "RIGID": {
                "type": "flag",
                "default": False,
                "description": "Rigid symmetry constraints",
            },
        },
    },

    # =========================================================
    # OUTPUT / WRITER to be replaced by G26
    # =========================================================
    "writer": {
        "panel": "Output",
        "description": "Preparation of input files for external codes",
        "keywords": {
            "GDV": {
                "type": "flag",
                "default": False,
                "description": "Prepare input for Gaussian 26",
            },
        },
    },

    # =========================================================
    # WORKFLOW / DEBUG
    # =========================================================
    "workflow": {
        "panel": "Workflow",
        "description": "Workflow control options",
        "keywords": {
            "DEBUG": {
                "type": "flag",
                "default": False,
                "description": "Enable debug mode",
            },
        },
    },
}
