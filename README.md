# Merlino 4.0

Merlino 4.0 starts from the closed Merlino 3.0 baseline and is the workspace for
the structural refactor. Merlino 3.0 remains frozen; new architecture work
happens here.

The first goal is to separate GUI, geometry/topology, Gaussian/GIC generation,
Fortran backends, DVR and data handling behind stable interfaces so each area
can evolve independently.

The planned scientific additions are:

- VPT2/VCI from canonical Merlino Hessian/QFF inputs, with Gaussian currently
  supported as an adapter and independent Python/Fortran77 kernels with
  Davidson diagonalization for large VCI spaces.
- Standard semiexperimental equilibrium geometries from Cartesian/GIC
  least-squares fits of isotopologue rotational data with QM vibrational
  corrections, avoiding fragile Z-matrix parameterizations.

See `doc/MERLINO4_REFACTOR_PLAN.md` before moving code.

## Main Entry Points

- `app.py`: GUI launcher.
- `python -m merlino`: Merlino4 workflow CLI for workspace initialization,
  GF/PED, VPT2/VCI, semiexperimental geometries and DVR command preparation.
- `python -m merlino_gui.app`: experimental Merlino4 workflow dashboard.
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
- `fortran/`: active Fortran backends. GICForge lives under
  `fortran/gicforge/`; the DVR kernel lives under `fortran/dvr/`; independent
  GF/VPT2/VCI/Davidson kernels live under `fortran/vpt2_vci/`; the
  semiexperimental geometry kernel lives under `fortran/semiexp/`.
- `bin/`: runnable binaries used by launchers, such as `gicforge.x`.
- `projects/`: local project/library data ignored by git.
- `working/`: runtime working directory ignored by git.
- `doc/`: architecture notes, freeze notes, reports, manuals, and current
  developer documentation.
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
under `fortran/gicforge/build/`; it updates `bin/gicforge.x` and the
compatibility alias `bin/prova.x`.

For the Fortran DVR backend:

```bash
cd fortran/dvr
./compile_MAC
```

The DVR build writes `bin/path_dvr.x`.

For the independent GF/VPT2/VCI source kernels:

```bash
cd fortran/vpt2_vci
./compile_check
```

This compiles `gf_core.f`, `vci_core.f` and `davidson_core.f` to object files.

For the semiexperimental geometry Fortran77 source kernel:

```bash
cd fortran/semiexp
./compile_check
```

## Merlino4 Runtime Contracts

New code should use the shared infrastructure in `merlino_core`:

- `MerlinoConfig` from `merlino.toml` for Gaussian/backend defaults.
- `WorkspaceLayout` for `inputs/`, `runs/`, `outputs/`, `reports/`, `cache/`
  and `logs/`.
- `RunManifest` for workflow reproducibility metadata and file checksums.
- Typed exceptions from `merlino_core.errors` for user-facing failures.

The CLI mirrors GUI-capable services without requiring Qt:

```bash
python -m merlino init my_project
python -m merlino gf --fchk calc.fchk --out gf_report.txt
python -m merlino vci --qff field.qff --max-quanta 3 --roots 6 --csv-dir csv
python -m merlino semiexp --xyz parent.xyz --observations isotopologues.csv --outdir semiexp_run --observable moments
python -m merlino gic --workdir working
python -m merlino gaussian-summary calc.log
python -m merlino backends
python -m merlino compare-backends
python -m merlino dvr-args --repo-root . --log scan.log --outdir out --figdir fig
```

See `doc/DEVELOPER_WORKFLOW.md` for the service/CLI/GUI/manifest contract used
by new workflows.

See `doc/SEMIEXPERIMENTAL_GEOMETRY.md` for the Merlino standard
semiexperimental geometry solver, including rationale, recommended defaults,
QM predicates and quality checks.
