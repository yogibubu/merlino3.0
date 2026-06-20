# Merlino4 GUI Design

## Direction

The Merlino4 GUI should be a workflow dashboard, not a collection of backend
buttons. It should collect inputs, display state and results, and call service
layers with documented contracts.

The legacy GUI under `gui/` remains available during migration. New GUI code
starts under `merlino_gui/`.

## Rules

- GUI classes must not contain scientific algorithms.
- GUI classes must not parse Gaussian logs directly.
- GUI classes must not call Fortran executables directly.
- GUI classes should call service interfaces and display normalized outputs.
- GUI workflow steps must ask for the backend family when alternatives exist
  (`python` or `fortran77`) and pass that choice to the service layer.
- Every long-running workflow should produce or consume a manifest JSON.

## Top-Level Workflows

- Molecule: structure input, isotopes, topology, rings and symmetry.
- GIC / Gaussian Input: GIC construction, GICForge and Gaussian input writing.
- DVR: Gaussian scan/path output to DVR levels, labels and reports.
- VPT2 / VCI: quartic force field to VPT2/VCI levels, with Davidson for large VCI.
- Semiexperimental Geometry: least-squares equilibrium geometry from
  isotopologue rotational constants and QM vibrational corrections.
- Jobs / Reports: manifests, logs, output files and reproducibility.

## Target Flow

```text
GUI -> service layer -> backend/parser/numerical code -> files + manifest -> GUI
```

The GUI should not implement backend-specific logic. It may expose a backend
choice to the user, but execution remains delegated to the workflow service,
which maps `python` or `fortran77` to the appropriate implementation.

## Migration Strategy

1. Add `merlino_gui` dashboard and workflow registry.
2. Keep `app.py` pointed at the legacy GUI until a workflow is fully migrated.
3. Move one workflow at a time behind a service interface.
4. Add result/report panels that read manifests instead of backend-specific
   temporary files.
5. Switch the main entry point only after GIC and DVR have service-backed GUI
   coverage.

## Experimental Entry Point

The new dashboard can be launched independently with:

```bash
python -m merlino_gui.app
```

It creates or reads `.merlino/project.json` under the selected working
directory. This project-state file is the first step toward making GUI state
explicit instead of implicit in individual windows.
