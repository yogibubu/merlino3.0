GUI — Merlino 3.0
Overview

The Merlino 3.0 GUI provides a minimal, controlled interface for preparing molecular inputs and managing updates to the central xyzin file.

Merlino4 adds a workflow dashboard on top of this legacy structure. The
dashboard is package-oriented: `merlino_gic` handles GIC definition/B matrices,
`merlino_gf` handles harmonic GF/PED, `merlino_vpt2_vci` handles anharmonic
VPT2/VCI, `merlino_semiexp` handles semiexperimental refinement and
`merlino_dvr` handles DVR. GUI controls build CLI/service calls and do not
duplicate scientific algorithms.

The GIC dashboard exposes symmetry adaptation as an explicit option at
definition time. If enabled, GIC names, irreps, point group and the
`symmetrized` flag are stored in `gic_definition.json`; B-matrix and GF/PED
steps receive this metadata from the schema rather than reassigning symmetry.

The semiexperimental dashboard exposes the same coordinate-model choice as the
CLI. `gic` is the default and uses GICForge-generated totally symmetric GICs.
`cartesian_symmetry` needs no Hessian and no GIC B matrix: it removes
translations/rotations from the Cartesian displacement space, symmetry-adapts
it, and fits only totally symmetric Cartesian directions. The final report
still writes primitive bond lengths, angles and dihedrals with propagated
errors.

The GUI is intentionally simple and follows a strict separation of responsibilities:

Input selection and preview → InputPanel

Commit confirmation and rollback → MainWindow

Visualization → ViewerPanel

Data exchange → exclusively via working/xyzin

Fundamental design rules

The GUI follows these non-negotiable rules:

The working directory is fixed:

<project root>/working


The only communication channel between modules is:

working/xyzin


GUI components never communicate directly with each other using internal data structures.

Any operation that modifies xyzin must:

be explicit

be confirmed by the user

support rollback

GUI structure
Main components
MainWindow
 ├── InputPanel
 └── ViewerPanel

Responsibilities
MainWindow

owns the application window

reacts to updates of xyzin

asks user confirmation after any update

performs rollback if the update is rejected

refreshes the viewer only after acceptance

MainWindow does not:

parse input files

generate structures

preview incomplete data

InputPanel

manages all user input

selects input type

handles preview vs commit

invokes readers

InputPanel does not:

ask confirmation

decide whether changes are accepted

manage rollback

ViewerPanel

visualizes the current molecular structure

reads data exclusively from working/xyzin

is refreshed only after confirmed updates

Input type selection

The input type is selected via horizontal radio buttons, not a dropdown menu.

Available input types:

[ SMILES ] [ XYZ ] [ Gaussian ] [ Molpro ] [ MRCC ]


This design:

avoids ambiguity

makes valid actions explicit

prevents invalid user gestures

SMILES input workflow
Preview vs commit

SMILES handling is explicitly split into two phases:

Preview

triggered on text change

non-destructive

does not modify xyzin

used only for live visualization

Commit

triggered only by pressing ENTER

writes to working/xyzin

invokes confirmation dialog via MainWindow

This guarantees:

no accidental updates

no repeated confirmation dialogs

full support for typing, paste, and drag-and-drop

File-based input workflow

File-based inputs include:

XYZ

Gaussian outputs

Molpro outputs

MRCC outputs

Workflow:

select input type

click Browse

select file

reader updates xyzin

confirmation dialog is shown

File inputs are always treated as commit operations.

Gaussian properties (separate file)
-----------------------------------
You can provide geometry from any supported source (SMILES/XYZ/MRCC/MOLPRO)
and a separate Gaussian LOG/OUT file for properties only. In this mode:

- The Gaussian file does NOT overwrite XYZ geometry.
- Only #VIBRATIONAL and #ROTATIONAL sections are updated.
- FCHK/FCH are not used in this pipeline (reserved for other workflows).

Decision logic (geometry vs harmonic vs anharmonic)
---------------------------------------------------
After input is committed, the GUI inspects #VIBRATIONAL:

- No frequencies → geometry-only workflow
- Frequencies without chi → harmonic workflow
- Frequencies + chi_cm1 → anharmonic workflow

Pipelines executed:
- Geometry-only: rot + thermo (trasl/rot only) + topology
- Harmonic: rot + thermo + DOS/Q(T) + topology
- Anharmonic: rot + thermo + DOS/Q(T) + topology

DOS/Q(T) settings
-----------------
The GUI provides a dedicated dialog to set:
- Emin and Emax (cm^-1) [always required]
- bin (cm^-1)
- vmax and ncap (for vibrational DOS)
- T (single temperature)
- emax_rot or jmax for rotational levels

Outputs (if VIBRATIONAL present):
- dos_vib.dat and vib_qt.dat
- dos_rovib.dat and rovib_qt.dat
For TS (negative frequency detected):
- n_vib_ts.dat (cumulative number of vibrational states)
- n_rovib_ts.dat (cumulative number of rovibrational states)

DOS/Q(T) caching
----------------
If `xyzin` and DOS settings are unchanged and the expected output files exist,
the GUI reuses cached DOS/Q(T) results and skips recomputation. Cache metadata
is stored in `working/dos_cache.json`.

Additional GUI artifacts:
- gui.log (timestamped events)
- summary.txt (input type, DOS settings, Q values, files generated)

Logging
-------
The GUI writes a single log file:

- working/gui.log

This log is shared across GUI components and records key events and workflow
steps (input updates, pipeline runs, BDPCS3 reports, warnings).

Workflow indicator
------------------
During pipeline execution the toolbar and input panel are disabled and the
status line shows a running message to prevent duplicate actions.

Refresh behavior
----------------
Viewer refresh is debounced to avoid repeated redraws during rapid updates.
The GUI also skips pipeline execution when `xyzin` is unchanged.

BDPCS3 report + Gaussian verification
-------------------------------------
The toolbar includes a "DPCS3 to BDPCS3 report" action. It:

- Prompts for BDPCS3 version: legacy or updated.
- Applies BDPCS3 bond corrections while keeping angles/dihedrals fixed.
- Writes a report with bond lengths, rotational constants (MHz), and coordinates.
- Generates companion files for verification.

Outputs in working/:
- bdpcs3.report
- bdpcs3.xyz
- bdpcs3.xyzin
- BDPCS3.gjf

After generation, the GUI shows quick actions to open/copy report and
Gaussian input paths.

The Gaussian input uses DPCS3 Cartesian coordinates and ModRedundant
bond targets (BDPCS3) so you can reproduce the same back-transform in
Gaussian (single-point + ModRedundant geometry update). The route line is:

#HF geom=modredundant output=pickett

Notes on Gaussian fchkin consistency
------------------------------------
If a working/fchkin file exists but has a different atom count than the
current working/xyzin, vibrations are skipped with a warning. This avoids
incorrect mass-weighting and prevents errors in downstream workflows.

Compact results view
--------------------
A small summary line below the status area shows Q_vib/Q_rovib and rotational
constants (A/B/C) when available.

Confirmation and rollback

Every successful update of xyzin triggers a confirmation dialog.

Options:

Yes → update accepted, viewer refreshed

No → previous xyzin restored from backup

This mechanism ensures:

no silent data loss

safe experimentation

reproducible workflows

Error handling

Preview errors are silently ignored

Commit errors propagate normally

GUI stability is always preferred over strict validation

The GUI is designed to never crash during preview.

## DVR Window

The main toolbar exposes a dedicated `DVR` action. It opens the Path DVR window
implemented in `advanced/dvr_window.py`.

The DVR window reads completed Gaussian scan/path logs and runs the vendored
backend `puckering_dvr/scripts/mw_path_dvr.py`. It does not generate Gaussian
inputs or scan paths; those remain part of the Merlino/`merlino_fit` Gaussian
preparation workflow.

Operational controls:

- `Use Latest Gaussian Log` selects the newest `.log`/`.out` file in the
  Merlino work directory.
- `Preflight / Preview` checks that the selected Gaussian file contains
  optimized geometries, energies or scan markers, and optional puckering GIC
  values.
- `Run Gaussian` starts Gaussian from the current `gauin.gjf`/`gauin`.
- `Run Gaussian Then DVR` runs Gaussian, selects the resulting log, previews it,
  and then starts the DVR analysis.
- `Refresh Results` previews the produced `*_summary.txt`, `*_levels.csv`,
  profile CSV and figures for the selected prefix.
- Each DVR run writes `*_dvr_run_manifest.json` with the exact command, Python
  executable, selected log, SHA256 checksum, output paths and Hamiltonian
  settings.

Default DVR settings:

- Gaussian log: `working/gauin.log`
- output directory: `working/puckering_dvr_outputs`
- figure directory: `working/puckering_dvr_figs`
- boundary: `periodic`
- solver: `fourier`
- rotational constants: enabled
- Cremer-Pople labeling: disabled unless explicitly requested

## GF / PED Window

The Advanced window exposes a dedicated `GF / PED` panel implemented in
`advanced/gf_window.py`.

The window reads a Cartesian Hessian from the current FCHK adapter, builds
Merlino non-redundant GICs or reuses a frozen GIC definition from
`gic-define`, evaluates the B matrix, transforms the Cartesian Hessian to the
GIC basis, optionally applies Pulay-style diagonal scaling with off-diagonal
geometric means, and reports frequencies, normal modes and PED.

GIC symmetrization is not done in this window. It must already be present in
the frozen definition produced by `gic-define`.

## VPT2 / VCI Window

The Advanced window exposes a separate `VPT2 / VCI` panel implemented in
`advanced/vpt2_vci_window.py`.

The window reads either a canonical indexed QFF text file or the current FCHK
adapter, applies mode/cutoff/pruning settings, and reports a VPT2/VCI energy
comparison plus dominant VCI contributions. This branch works in Cartesian
normal modes and does not use GICs.

Both windows support `Export Report`, `Export CSVs`, `Save Preset` and
`Load Preset`. The GUI does not contain scientific algorithms; it only
validates paths and dispatches report export/preset actions. Each successful
GF or VPT2/VCI GUI run writes a `merlino.run.v1` manifest in the current work
directory.

GUI testing policy

Core GUI widgets are covered by focused Qt tests under `gui/tests`. Full
interactive workflows are still checked manually with:

```bash
python app.py
```

focus automated tests on scientific correctness

Scope and limitations

At this stage, the GUI:

does not manage projects

does not store session history

does not perform calculations

does not integrate external editors unless explicitly added

All such features are considered out of scope for Merlino 3.0.

Final note

The Merlino 3.0 GUI is designed to be:

predictable

explicit

minimal

scientifically safe

Every GUI decision prioritizes data integrity over convenience.
