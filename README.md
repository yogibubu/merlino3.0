# Merlino 4.0

Merlino 4.0 starts from the closed Merlino 3.0 baseline and is the workspace for
the structural refactor. Merlino 3.0 remains frozen; new architecture work
happens here.

The first goal is to separate GUI, geometry/topology, Gaussian/GIC generation,
Fortran backends, DVR and data handling behind stable interfaces so each area
can evolve independently.

See `doc/MERLINO4_REFACTOR_PLAN.md` before moving code.

## Main Entry Points

- `app.py`: GUI launcher.
- `manager.py`, `cli_modules.py`, `cli_gaussian.py`: command-line entry points
  retained for compatibility.
- `freeze_check.sh`: regression check before a freeze or commit.
- `freeze_pack.sh`: creates a source archive after validation.

## Functional Folders

- `gui/`: main PySide6 GUI, input handling, viewers, and GUI tests.
- `advanced/`: advanced calculation windows and launchers for Gaussian, GICForge,
  DVR, VPT2, MSR, and `merlino_fit`.
- `geometry/`: rotational, vibrational, thermodynamic, rovibrational, symmetry,
  and `xyzin` geometry utilities.
- `merlino_fit/`: topology, synthons, BDPCS3/survibfit workflows, Gaussian GIC
  generation, and related tests.
- `puckering_dvr/`: vendored DVR backend. It consumes completed Gaussian outputs
  and does not generate Gaussian paths.
- `fortran/`: Fortran backends and legacy numerical programs. Active GICForge
  source is under `fortran/gicforge/`.
- `bin/`: runnable binaries used by launchers, such as `gicforge.x`.
- `projects/`: local project/library data ignored by git.
- `working/`: runtime working directory ignored by git.
- `doc/`: architecture notes, freeze notes, reports, manuals, and legacy root
  notes.
- `scripts/`: smoke scripts and maintenance utilities.

## Validation

Run the inherited validation suite:

```bash
PYTHON=python ./freeze_check.sh
```

For the Fortran GICForge backend:

```bash
cd fortran/gicforge
./compile_MAC
```

The Fortran compile script uses legacy-compatible flags and writes build logs
under `fortran/gicforge/build/`; it updates `fortran/gicforge/gicforge`,
`bin/gicforge.x`, and the compatibility alias `bin/prova.x`.

For the Fortran DVR backend:

```bash
cd fortran/dvr
./compile_MAC
```

The DVR build writes `bin/path_dvr.x`.
