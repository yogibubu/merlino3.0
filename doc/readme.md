# Merlino 3.0 Developer Map

Merlino is organized as a set of independent scientific packages that exchange
state through a common `xyzin` molecular container.  New workflows should avoid
private side channels: geometry, topology, rotational, vibrational and
isotopologue data belong in explicit `xyzin` sections or in run manifests.

## Runtime Packages

- `merlino_core`: shared configuration, manifests, `xyzin` section handling,
  isotopologue records and numerical helpers.
- `merlino_gic`: GIC definition, symmetry assignment and reusable B-matrix
  evaluation.
- `merlino_semiexp`: SEfit preprocessing, constraints, least-squares solver,
  diagnostics and reports.
- `merlino_gf`: Cartesian Hessian reading, GIC/GF transformation, frequencies,
  normal modes and potential-energy distributions.
- `merlino_dvr`: DVR execution helpers for completed Gaussian scans/paths.
- `merlino_vpt2_vci`: normal-mode quartic force-field VPT2/VCI utilities.
- `merlino_gui`: new workflow dashboard and manifest browser.
- `gui` and `advanced`: compatibility GUI panels still used during migration.
- `fortran`: active Fortran77-compatible numerical backends.

## `xyzin` Contract

`xyzin` starts with a standard XYZ block and then contains uppercase sections
such as `#BASIC`, `#ROTATIONAL`, `#VIBRATIONAL`, `#TOPOLOGY` and
`#ISOTOPOLOGUES`.  A module that regenerates a section must preserve unrelated
sections.  The current section format is documented in
`doc/XYZIN_FORMAT.md`.

SEfit always runs from `xyzin`.  External CSV/JSON/TOML/MSR inputs are
preprocessing sources: they are first materialized into the common container and
then reread from there.

## Runtime And Generated Files

- `working/` and `tmp/` are runtime-only and ignored by git.
- `working/xyzin` is a local workspace file, not a versioned source artifact.
- Build intermediates belong under ignored build directories.
- Scientific examples and reproducibility inputs belong under `examples/` or
  `benchmarks/`, not under `working/`.

## Current Documentation

- `doc/REPOSITORY_LAYOUT.md`: repository map and build policy.
- `doc/PACKAGE_ARCHITECTURE.md`: package boundaries and data flow.
- `doc/XYZIN_FORMAT.md`: canonical `xyzin` and isotopologue section format.
- `doc/GICFORGE.md`: current GICForge behavior.
- `doc/GICFORGE_HBOND_FRAGMENT_ROADMAP.md`: planned H-bond and fragment
  coordinate extensions.
- `doc/SEMIEXPERIMENTAL_FILE_FORMATS.md`: SEfit external compatibility formats.
- `doc/SEMIEXPERIMENTAL_GEOMETRY.md`: SEfit methodology and outputs.
- `doc/testing.md`: regression and verification workflow.

Historical cleanup, freeze and triage notes are kept under `doc/archive/`.
