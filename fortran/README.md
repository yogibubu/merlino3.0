# Fortran Source Tree

This tree contains only active Fortran backends used by Merlino4. Generated
executables do not belong in source directories.

Active backends:

- `gicforge/`: Cartesian XYZ to GIC/Gaussian input backend.
- `dvr/`: Fortran77 DVR numerical kernel.
- `harmonic_internal/`: Wilson GF harmonic internal-coordinate source backend.

Build products:

- `fortran/gicforge/build/gicforge`
- `fortran/dvr/build/path_dvr`
- `fortran/harmonic_internal/build/gf.o`
- runtime copies under `bin/`

Historical standalone Fortran material is kept in the frozen Merlino3.0 tree,
not duplicated here.

The root of this tree should stay descriptive only. Put active code in a
backend subdirectory.

Keep generated binaries, logs, object files and temporary compiler outputs out
of this tree unless they are explicitly promoted to runtime binaries under
`bin/`.
