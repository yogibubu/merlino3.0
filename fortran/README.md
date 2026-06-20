# Fortran Source Tree

This tree contains source code and small reference inputs only. Generated
executables do not belong in source directories.

Active backends:

- `gicforge/`: Cartesian XYZ to GIC/Gaussian input backend.
- `dvr/`: Fortran77 DVR numerical kernel.

Build products:

- `fortran/gicforge/build/gicforge`
- `fortran/dvr/build/path_dvr`
- runtime copies under `bin/`

Historical or reference material:

- `legacy/`: old vibrational, VCI/DVR and rate utilities.
- `gnic/`: standalone GNIC development experiments.
- `qcent/`, `quadrupolari/`, `symmetry/`, `volt/`: reference utilities and data.

Keep generated binaries, logs, object files and temporary compiler outputs out
of this tree unless they are explicitly promoted to runtime binaries under
`bin/`.
