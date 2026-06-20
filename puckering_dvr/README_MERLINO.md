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
- run Gaussian
- run `Puckering DVR – Gaussian scan analysis`

The DVR panel reads the Gaussian log, extracts the optimized scan structures,
builds the mass-weighted path coordinate and writes levels, profiles and plots
under the selected output directories.

## Command-Line Equivalent

```bash
python puckering_dvr/scripts/mw_path_dvr.py \
  --gaussian-log working/gauin.log \
  --log-selection last-per-link \
  --boundary periodic \
  --solver fourier \
  --compute-rotconst \
  --label-cremer-pople \
  --outdir working/puckering_dvr_outputs \
  --figdir working/puckering_dvr_figs \
  --prefix puckering_dvr
```

For the Merlino ring-GIC path, Gaussian uses inactive `RPck....` coordinates
and active `QPck....`/`PhiP....` coordinates. `PhiP` is interpreted as a phase
angle and `QPck` as the torsional puckering amplitude used by the current GIC
parametrization. The DVR post-processing can label the path with
Cremer-Pople-like quantities when requested.
