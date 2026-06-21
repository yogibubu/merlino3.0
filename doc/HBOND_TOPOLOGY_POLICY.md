# Hydrogen-Bond Topology Policy

Merlino keeps covalent topology and non-covalent interactions separate.
This prevents directional hydrogen bonds from creating pseudo-rings in GIC
construction before the redundancy treatment for non-covalent cycles has been
defined.

## Graphs

- Covalent graph: used by GICForge for primitive coordinates, ring detection,
  redundancy removal, symmetry adaptation and Gaussian GIC generation.
- Interaction graph: used for non-covalent targets such as X-H...Y hydrogen
  bonds in geometry correction, constraints and reporting.

Hydrogen bonds must not be inserted into the covalent graph by default.

## BDPCS3 Geometry Correction

BDPCS3 may add an H...Y distance target for a directional X-H...Y contact.
The target is a non-covalent correction coordinate, not a GIC primitive:

```text
kind = hbond_target
topology = noncovalent
participates_in_gic = false
participates_in_ring_detection = false
```

The internal-to-Cartesian back-transform may include this target with its own
metric weight. Other primitive coordinates remain covalent-topology primitives.

## Shared Parameters

The source of truth is:

```text
merlino_core/parameters/bdpcs3_hbond.toml
```

Python reads this file through `merlino_core.parameters.bdpcs3`. Fortran77 uses
the synchronized include:

```text
merlino_core/parameters/fortran/bdpcs3_hbond_params.inc
```

The test suite checks that the TOML and Fortran include contain the same
numeric values.

Regenerate the Fortran include after changing the TOML with:

```bash
python -m merlino_core.parameters.generate_fortran_includes
```

Current defaults:

- X-H-Y angular threshold: 150 deg.
- H...Y distance damping: error function centered at 3.0 Angstrom, width
  0.15 Angstrom.
- H-bond search cutoff: 3.8 Angstrom.
- O-H...O correction: -0.055 Angstrom.
- N-H...N correction: -0.055 Angstrom.
- Mixed N/O donor-acceptor corrections: average of the N-H...N and O-H...O
  values.

## GICForge

`FindHBnd` detects hydrogen bonds without modifying `NBond/IBond`.
`MkHBnd` remains as a legacy wrapper capable of adding contacts to topology,
but the standard GICForge workflow calls `FindHBnd` only. Therefore H-bonds are
reported as non-covalent targets and do not generate pseudo-cycles.

If non-covalent GICs are added later, they must be implemented as a separate
coordinate block. Redundancy removal must occur inside that block first, and a
non-covalent coordinate must not replace a covalent stretching, bending or
torsional coordinate.
