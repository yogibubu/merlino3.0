# Puckering DVR in Merlino 3.0

This directory is the operational copy of the puckering DVR workflow imported
from:

`/Users/vincenzobarone/Desktop/puckering_dvr_github`

The Merlino copy intentionally includes code, documentation and runnable
examples, but excludes the source `.git` directory, the legacy zip package and
the bibliography PDFs. Those files are not needed at runtime.

## Merlino Flow

1. Start the environment:

```bash
source ~/.bashrc
merlino-set
merlino-run-check
```

2. Open the GUI:

```bash
merlino-run
```

3. In Advanced Calculations:

- generate the Gaussian input from the current `xyzin`

4. In the toolbar `DVR` window:

- run Gaussian from `gauin.gjf`, or select an already completed Gaussian log
- run preflight/preview
- run DVR
- inspect the summary, levels CSV, profile CSV and figures

The DVR panel reads any Gaussian log containing optimized scan/path structures.
The scanned coordinate does not have to be puckering: the Hamiltonian coordinate
is always the mass-weighted Cartesian path length built from consecutive
optimized geometries. Ring/Cremer-Pople labels are optional post-processing
columns.

Each GUI DVR run writes `*_dvr_run_manifest.json` beside the DVR output files.
The manifest records the exact command line, Python executable, selected
Gaussian log, SHA256 checksum, output directories, prefix, boundary condition,
solver and labeling flags.

## Command-Line Equivalent

```bash
python puckering_dvr/scripts/mw_path_dvr.py \
  --gaussian-log working/gauin.log \
  --log-selection last-per-link \
  --boundary periodic \
  --solver fourier \
  --compute-rotconst \
  --outdir working/puckering_dvr_outputs \
  --figdir working/puckering_dvr_figs \
  --prefix puckering_dvr
```

Add `--label-cremer-pople --ring 1,2,3,4,5` only when the Gaussian scan is a
ring-puckering scan and Cremer-Pople labels are wanted.

For the Merlino ring-GIC path, Gaussian uses inactive `RPck....` coordinates
and active `QPck....`/`PhiP....` coordinates. `PhiP` is interpreted as a phase
angle and `QPck` as the torsional puckering amplitude used by the current GIC
parametrization. When those GIC values are present in the Gaussian log, the DVR
profile writes the raw `gic_QPck....`/`gic_PhiP....` values and, if
Cremer-Pople labeling is enabled, a fitted bridge from Gaussian puckering
components to Cartesian Cremer-Pople components.
