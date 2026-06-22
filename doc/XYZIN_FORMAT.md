# Merlino `xyzin` Container

`working/xyzin` is the common Merlino molecular container used by the older
synthons, rotational and vibrational modules and by the newer GIC/GF/SEfit
workflow.

The file starts with a standard XYZ block. Additional data are appended as
uppercase sections. A module that regenerates one kind of data must replace only
its own section and preserve all unrelated sections.

Typical sections are:

- `#BASIC`: charge, spin multiplicity, point group and thermodynamic defaults.
- `#SMILES`: source SMILES, when the input came from a SMILES reader.
- `#ROTATIONAL`: parent rotational constants, DeltaVib bridge values and
  rotational metadata.
- `#VIBRATIONAL`: harmonic/anharmonic vibrational data.
- `#THERMO`: thermodynamic functions.
- `#GAUSSIAN_TOPOLOGY`: optional Gaussian-derived charges and bond orders.
- `#TOPOLOGY`: topology and symmetry data.
- `#ISOTOPOLOGUES`: isotopic definitions, optionally enriched with rotational
  observations and corrections.

## `#ISOTOPOLOGUES`

The isotopologue section has schema `merlino.xyzin.isotopologues.v1`. It is a
general Merlino section, not a SEfit-only format. Each record must define the
isotopologue label and isotopic substitution. Rotational constants,
vibrational corrections, electronic corrections and standard deviations are
optional enrichment lines. Modules that know only isotope substitutions write
definition-only records; modules that know the spectroscopic data add the
corresponding optional lines in the same block.

Example:

```text
#ISOTOPOLOGUES
SCHEMA merlino.xyzin.isotopologues.v1
UNITS ROTATIONAL=MHz DELTAVIB=MHz DELTAEL=MHz SIGMA=MHz
INDEXING ATOMS=ONE_BASED
BEGIN parent
DEFINITION parent
ROTATIONAL_MHZ A=1000 B=800 C=600
DELTAVIB_MHZ A=1 B=2 C=3 SOURCE=gaussian CONVENTION=subtract
DELTAEL_MHZ A=0 B=0 C=0 SOURCE=unspecified CONVENTION=subtract
SIGMA_MHZ A=0.01 B=0.02 C=0.03
END
BEGIN D2
DEFINITION 2:2
ROTATIONAL_MHZ A=990 B=790 C=590
DELTAVIB_MHZ A=0.5 B=0.6 C=0.7 SOURCE=gaussian CONVENTION=subtract
DELTAEL_MHZ A=0 B=0 C=0 SOURCE=unspecified CONVENTION=subtract
END
BEGIN D3_definition_only
DEFINITION 3:2
END
```

`DEFINITION` uses one-based atom indexes and compact substitutions:
`atom_index:mass_number`, separated by semicolons for multiple substitutions,
for example `2:2;5:13`. `parent` means no substitution.

SEfit starts from `xyzin`. External TOML/JSON/CSV/MSR/job-table data are
preprocessing sources only: they are first materialized into `#ISOTOPOLOGUES`
and then the fit rereads geometry and observations from the updated `xyzin`.
The CLI option `--xyzin` selects the canonical container explicitly; the GUI
always passes its project `xyzin`. If a `xyzin` file does not exist, the
preprocessing step creates one from the supplied Cartesian/MSR/job geometry and
appends the isotopologue section. SEfit can run only when the relevant records
contain complete `ROTATIONAL_MHZ` lines.

Validation is centralized in `merlino_core.isotopologues`.  The shared
validator checks schema version, duplicate labels, duplicate definitions,
one-based atom indexes, atom-count compatibility when a geometry is available,
finite positive rotational constants, finite corrections, positive
experimental standard deviations and supported correction conventions.  Modules
that need only isotope definitions can accept definition-only records; SEfit
calls the same validator with `require_rotational=True`.
