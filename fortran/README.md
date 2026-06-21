# Fortran Source Tree

This tree contains only active Fortran backends used by Merlino4. Generated
executables do not belong in source directories.

Active backends:

- `gicforge/`: Cartesian XYZ to GIC/Gaussian input backend.
- `bdpcs3/`: Fortran77 mirror of the Python topology/synthon BDPCS3 backend.
- `dvr/`: Fortran77 DVR numerical kernel.
- `vpt2_vci/`: independent GF, VCI and Davidson source kernels.

Build products:

- `fortran/gicforge/build/gicforge`
- `fortran/bdpcs3/build/bdpcs3`
- `fortran/dvr/build/path_dvr`
- `fortran/vpt2_vci/build/*.o`
- runtime copies under `bin/`

`gicforge/mkprim.f` also provides weighted metric builders for
internal-to-Cartesian back-transform (`MakeGW`, `MkGm1BW`). The weights are
least-squares metric weights on internal coordinates; primitive generation is
unchanged.

Hydrogen bonds are detected as non-covalent correction targets, not inserted
into the covalent topology used by GICForge. The shared policy is documented in
`doc/HBOND_TOPOLOGY_POLICY.md`, and the shared numeric parameters live under
`merlino_core/parameters/`.

Historical standalone Fortran material is kept in the frozen Merlino3.0 tree,
not duplicated here.

The root of this tree should stay descriptive only. Put active code in a
backend subdirectory.

Keep generated binaries, logs, object files and temporary compiler outputs out
of this tree unless they are explicitly promoted to runtime binaries under
`bin/`.
