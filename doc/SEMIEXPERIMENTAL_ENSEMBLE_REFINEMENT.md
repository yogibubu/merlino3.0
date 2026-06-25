# Multi-molecule class-constrained semiexperimental refinement

This note defines the first MORPHEUS ensemble-refinement layer.  It is intended
for homologous or chemically related molecules where isotopic substitution is
too sparse to support a full independent semiexperimental structure for each
system, but where the systematic error of the computational reference geometry
is expected to be transferable within chemically defined coordinate classes.

## Model

For molecule `m` and generated coordinate `k`, the ensemble model is

```text
q_SE(m,k) = q_QC(m,k) + Delta[class(m,k)]
```

The absolute coordinate values are not forced to be equal across molecules.
Only the correction relative to the computational reference is shared.  This is
the physically defensible statement for large related molecules: two C-C bonds
in different steroids need not be equal, but the transferable QC bias of a
specified C-C class can be refined from all available rotational data at once.

Separate classes can be used for chemically distinct families, for example
`CC_short`, `CC_long`, `CO_single`, `CO_carbonyl`, `CH_methyl`, or
angular/torsional families.  A class must not mix coordinate types; the
implementation rejects class definitions that match, for example, both
stretches and bends.

## Current implementation

The public API is:

```python
from merlino_semiexp import (
    EnsembleClassCorrection,
    EnsembleMolecule,
    fit_ensemble_class_corrections,
)

result = fit_ensemble_class_corrections(
    (
        EnsembleMolecule("testosterone", testosterone_request),
        EnsembleMolecule("androsterone", androsterone_request),
    ),
    (
        EnsembleClassCorrection("CC_skeleton", ("R(", "C,C")),
        EnsembleClassCorrection("CO_carbonyl", ("C=O",)),
    ),
)
```

Each molecule remains an ordinary `SemiexperimentalFitRequest` with its own
geometry, topology, symmetry, isotopologues, vibrational corrections, component
selection, and coordinate model.  The ensemble layer builds the generated GIC
model for each molecule and constructs a single weighted linear least-squares
problem for the shared class corrections.

The same model can be run from a TOML job:

```bash
python -m merlino_core.cli semiexp-ensemble \
  --job examples/semiexp/anhydrides_ensemble/anhydrides_ensemble.mse-ensemble.toml \
  --outdir working/semiexp/anhydrides_ensemble
```

The prior policy can be benchmarked systematically with:

```bash
python -m merlino_core.cli semiexp-ensemble-compare \
  --job examples/semiexp/anhydrides_ensemble/anhydrides_ensemble.mse-ensemble.toml \
  --outdir working/semiexp/anhydrides_prior_comparison
```

This writes three variants: `no_prior`, `soft_prior`, and `hard_constraint`.
The hard-constraint variant removes angular class corrections from the global
variable set, which is equivalent to fixing those corrections to zero.
The reported condition number is computed after column scaling of the weighted
design matrix; fitted corrections, covariance and standard deviations are then
back-transformed to the original coordinate units.

The manuscript artifacts for the ensemble-refinement JPCL skeleton can be
regenerated with a single command:

```bash
python -m merlino_core.cli semiexp-ensemble-paper \
  --job examples/semiexp/anhydrides_ensemble/anhydrides_ensemble.mse-ensemble.toml \
  --paper-dir doc/papers/ensemble_jpcl \
  --outdir working/semiexp/anhydrides_full_analysis
```

This command refreshes the comparison CSV files, prior-strength scan,
leave-one-molecule-out diagnostics, and the LaTeX fragments used by the paper.
All ensemble CLI entry points also write a canonical `run_manifest.json` with
input/output checksums.  This is separate from the scientific
`ensemble_manifest.json`: the former freezes a workflow execution, while the
latter describes the fitted model.

## Stable ensemble job schema

An ensemble job uses the schema identifier
`merlino.semiexp.ensemble.v1`.  The stable top-level sections are:

- `[fit]`: numerical finite-difference and rank threshold controls.
- `[acceptance]`: rank, conditioning, support and correlation policy.
- `[[molecules]]`: ordinary single-molecule SE jobs or legacy MSR inputs.
- `[[classes]]`: transferable coordinate-class definitions and optional
  Gaussian soft priors.

The default acceptance policy is deliberately conservative and should be
written explicitly in production examples:

```toml
[acceptance]
require_full_rank = true
max_condition_number = 1.0e8
min_residual_degrees_of_freedom = 1
min_molecule_support = 2
high_correlation_review_threshold = 0.98
high_correlation_reject_threshold = 0.9999
```

The result status is `accepted`, `review`, or `rejected`.  `accepted` means the
model passes the rank, conditioning, molecule-support and correlation checks.
`review` means the model is numerically usable but contains strongly
correlated class corrections that should be interpreted explicitly.  `rejected`
means the correction model should not be used as a production model without
changing classes, priors or atom typing.  The status is written to the text
report, CSV summaries and `ensemble_manifest.json`.

The acceptance policy is part of the scientific manifest.  A production fit is
therefore defined by the input job, the generated class projections, the
weighted design matrix diagnostics, and the exact policy that accepted or
rejected the result.  The minimum policy for a transferable model is:

- full rank unless an explicitly reduced hard-constraint model is being tested;
- finite scaled condition number below the declared threshold;
- at least one residual degree of freedom;
- support from the requested number of molecules for every transferable class;
- explicit review or rejection of highly correlated class corrections.

The ensemble job points to ordinary single-molecule SE job files.  Those job
files may keep the Cartesian reference geometry and isotopologue table in
separate files:

```toml
schema = "merlino.semiexp.job.v1"

[files]
geometry = "parent.xyz"
observations = "isotopologues.toml"
```

This keeps literature-data cases portable and avoids duplicating large
coordinate blocks in the fit definition.

Class selectors can be substring based, coordinate-kind based, atom-symbol
based, primitive-value based, synthon-based, or a combination.  For
symmetrized GICs, a class correction is projected through the primitive
coefficients of the active coordinate rather than treated as an all-or-nothing
label match.  Thus a shared `CC_short` correction can contribute to an A1
coordinate that contains several C-C primitive stretches, but only those whose
reference distance falls inside the declared value window.

Atom-symbol selectors are type aware:

- stretches use the unordered atom pair;
- bends preserve the central atom and canonicalize only the two terminal atoms;
- torsions use the unordered central-atom pair, i.e. the chemistry of the
  central bond;
- out-of-plane coordinates use the central atom of `U(center,a,b,c)`.

Soft priors are declared directly on classes:

```toml
[[classes]]
name = "CCO_bend"
kind = "bend"
atoms = ["C", "C", "O"]
prior_value = 0.0
prior_sigma = 1.0e-3
```

They are added as Gaussian observations on the class correction
`Delta_class = prior_value +/- prior_sigma`.  They regularize weakly supported
shared corrections while leaving the spectroscopic residuals and per-molecule
diagnostics separate.

Stretch selectors can be split by reference primitive values:

```toml
[[classes]]
name = "CO_carbonyl"
kind = "stretch"
atoms = ["C", "O"]
value_max = 1.28

[[classes]]
name = "CO_single"
kind = "stretch"
atoms = ["C", "O"]
value_min = 1.28
```

The preferred synthon selector is continuous.  Atomic synthons are condensed
into the effective atomic number `Zeff`; two atoms are assigned to the same
type when their `Zeff` values differ by less than a declared threshold.  For a
primitive coordinate, the comparison uses the chemically relevant ordered
signature: unordered pair for stretches, central atom retained for bends,
central bond for torsions, and central atom for out-of-plane coordinates.

```toml
[[classes]]
name = "Ccarb_Ocarbonyl"
kind = "stretch"
atoms = ["C", "O"]
synthon_zeff = [5.91, 7.865]
synthon_threshold = 0.035
```

Discrete `synthon_signatures` remain available for reproducibility checks, but
they are less expressive than thresholded `Zeff` atom types.

The first implementation is deliberately linearized around the computational
reference geometries.  It returns:

- global class corrections and standard errors;
- covariance and correlation between global class corrections;
- an explicit acceptance status and the reasons/review items behind it;
- weighted RMS before and after the shared correction;
- per-molecule row counts, matched-coordinate counts, rank contribution, and
  residual reduction.
- leave-one-molecule-out transfer scores, defined as the predicted held-out
  WRMS divided by the uncorrected held-out WRMS.

The output writer creates:

- `ensemble_class_corrections.txt`
- `ensemble_class_corrections.csv`
- `ensemble_class_report.csv`
- `ensemble_molecule_blocks.csv`
- `ensemble_manifest.json`
- `ensemble_covariance.csv`
- `ensemble_correlation.csv`
- `run_manifest.json` when the calculation is launched from the CLI

The implementation fails early if any declared class matches no active
coordinate.  Rank and condition number are reported explicitly.  For the three
anhydrides, the robust demonstration model uses shared C-C, C-O and C-H stretch
corrections plus softly regularized bend classes.  Without priors the broad
bend classes are rank-consistent but poorly conditioned, which is precisely the
case soft priors are meant to handle.

## Regression set

The current automated regression set is intentionally small and scientific:

- anhydrides ensemble: no-prior, soft-prior and hard-constraint variants;
- glycine I/II ensemble: continuous synthon `Zeff` atom typing and threshold
  scan;
- primitive signature canonicalization: torsions are matched by central bond,
  and out-of-plane coordinates by central atom;
- GUI job-writer round trip: fine selectors, soft priors and absolute paths are
  preserved;
- CLI run manifest: anhydrides ensemble writes both scientific and workflow
  manifests.

Run the focused checks with:

```bash
pytest merlino_semiexp/tests/test_ensemble.py gui/tests/test_ensemble_window_helpers.py
```

## Scientific use

This layer is meant for parent-only or isotope-poor large molecules, especially
homologous series such as steroid hormones.  It replaces isotopic richness with
chemical-class richness across molecules.  A future short paper could focus on
this point separately from the main MORPHEUS manuscript:

1. define the transferable correction model;
2. validate it on small molecules where isotope-rich benchmarks exist;
3. demonstrate the scale case on a hormone series;
4. compare broad classes (`CC`) against chemically separated classes (`CC`,
   `CO_single`, `CO_carbonyl`, etc.);
5. quantify robustness from adding molecules to the same global fit.

## Next steps

- Add a nonlinear polishing stage that applies the class corrections to each
  molecule and optionally relinearizes.
- Add larger homologous series as soon as complete experimental constants and
  vibrational corrections are available.
