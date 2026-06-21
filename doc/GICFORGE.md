# GICForge

GICForge is the active Fortran77 geometry backend formerly called `prova`.

## Scope

GICForge is intentionally narrow:

- read Cartesian molecular input only from `xyzin` in the Merlino working
  directory
- determine the molecular point group through the linked Fortran `symm.f`
  symmetry engine
- build topology, primitive stretch coordinates and redundant non-stretch GIC
  candidates
- reduce non-stretch candidates to non-redundant GICs by coordinate-type
  blocks, so bends, linear bends, torsions and out-of-plane coordinates are
  never mixed during the final rank pruning
- optionally symmetrize comparable GICs inside each coordinate-type block only
- canonicalize ring atom numbering with the same Prelog-first convention used
  by Python
- optionally build/write the B matrix used for geometry-parameter changes
- write a readable geometry/GIC report
- write Gaussian input with inactive `RPck....` puckering components and active
  `QPck....`/`PhiP....` coordinates

Everything else belongs to Python: GUI orchestration, RDKit/SMILES, project
management, DVR, Cremer-Pople post-processing, regression comparison, freeze
checks and user-facing workflow logic.

## Utility Architecture

The Merlino4 GIC layer is split into reusable utilities.

1. `gic-define` receives Cartesian coordinates and atomic symbols/numbers,
   runs GICForge once, optionally performs the deterministic symmetry post-check, and
   writes a frozen `merlino.gic.definition.v1` JSON schema.  The schema stores
   the primitive coordinates, the GIC coefficient matrix, labels, irreducible
   representations, point group, a `symmetrized` flag and the
   Gaussian-readable GIC block. Use `--no-symmetry` when raw non-redundant GICs
   must be frozen without symmetry adaptation.
2. `gic-bmatrix` receives a frozen GIC schema plus a current Cartesian
   geometry and atomic symbols/numbers, and evaluates GIC values and the Wilson
   B matrix.  It also propagates the frozen GIC names, irreps, point group and
   `symmetrized` flag to downstream programs or to an optional metadata CSV.
   It does not perform topology perception, redundancy removal or symmetry
   assignment.
3. `gic-gf` receives a frozen GIC schema plus a Cartesian Hessian adapter
   currently implemented for Gaussian FCHK. It evaluates the same frozen GICs
   on the Hessian geometry, transforms the Cartesian Hessian to the GIC basis,
   optionally applies Pulay-style internal-Hessian scaling, solves Wilson GF,
   and writes frequencies, normal modes, G matrix, internal force constants and
   PED tables.

The second utility is also exposed as a Python library function,
`merlino_gic.evaluate_gic_definition`, so other programs can build B matrices
without launching GICForge.  The supplied atomic numbers only need to be
compatible in atom count with the frozen schema; the B matrix itself is defined
by the primitive coordinate expressions and the current Cartesian coordinates.
The returned evaluation object includes `labels`, `names`, `irreps`,
`point_group` and `symmetrized` metadata in addition to values and B matrices.
For fit/diagnostic use, `gic-bmatrix` reports derivatives in the units of the
supplied geometry file. For Hessian transformation, `gic-gf` evaluates B with
Cartesian coordinates in bohr, matching the canonical Hessian units
Eh/bohr^2.

In semiexperimental refinements, `SEfit` calls the definition utility only at
the beginning of a fresh fit.  All ordinary iterations reuse the same frozen
GIC schema and rebuild only the B projector when necessary.  A restart is the
explicit boundary at which the GIC schema may be regenerated.

## Symmetrization Ownership

GIC symmetrization belongs exclusively to `gic-define`. The sequence is:

1. GICForge constructs primitive and non-redundant coordinates.
2. GICForge and the deterministic Python post-check assign the point group,
   irreducible representations and symmetry-adapted final GIC labels.
3. The resulting `merlino.gic.definition.v1` schema is frozen.

Downstream programs must not re-symmetrize. `gic-bmatrix`, `gic-gf`,
semiexperimental refinement and Gaussian input writers only evaluate, filter
or order the frozen coordinates. If a different symmetry tolerance, geometry
or coordinate policy is required, the correct operation is a new `gic-define`
run or an explicit restart that creates a new schema.

Pulay scaling files accepted by `gic-gf` are line-oriented text or CSV files:

```text
# selector factor
default 1.000
GIC003 0.980
A1Str0001,0.995
5=1.020
```

Selectors may be `default`/`all`, one-based indices, `GICnnn` labels, exact
GICForge names, exact labels, or unambiguous substrings. If `s_i` and `s_j`
are the diagonal factors, the internal Hessian element is scaled as
`F_ij <- F_ij sqrt(s_i s_j)`.

The Fortran GICForge build no longer contains a SMILES reader and no longer
accepts legacy FITPOT/VCI/DVR, MSR/isotope or rate keywords. Those workflows
must be driven from Python or from explicitly archived legacy code.

## Input Contract

GICForge no longer selects among geometry readers. The only molecular geometry
input is the XYZ-format file named `xyzin` in the run directory. The `provin`
file is still used for calculation keywords, title, charge and multiplicity,
but it must not contain Cartesian coordinates, Z-matrices or FCHK geometry
instructions.

Unsupported legacy geometry paths:

- FCHK geometry input
- Z-matrix input from `provin`
- Cartesian coordinates embedded in `provin`
- SMILES input or SMILES-derived coordinate generation
- `BLDZ`/`WRTZ` Z-matrix build/write workflow
- `FITPOT`/`VCI1`/`DVR1`, `SEMIEX`/`ISOTOP`/`DVIBROT`, `VOLT`, `COLL`/`GORIN`
  and `SPINFOR` in the GICForge driver

`RDXYZ` is therefore implicit and is no longer a user-facing keyword.

## Build

```bash
cd fortran/gicforge
./compile_MAC
```

The build creates:

- `fortran/gicforge/build/gicforge`
- `bin/gicforge.x`
- `bin/prova.x` as a compatibility alias

The source remains fixed-form Fortran77/legacy-compatible Fortran.

## Molecular Symmetry

GICForge links `fortran/gicforge/symm.f` directly in the normal build. After
Cartesian input and mass data are available, `coord.f` centers the molecule at
the center of mass, projects it on the principal inertia axes, and calls
`DETERMINE_POINT_GROUP`.

The readable output reports:

- `Point Group from symm.f`
- symmetry quality: `STRICT`, `QUASI`, or `BROKEN`
- maximum atom-matching deviation in Angstrom

The second line of `xyzin` may still contain an optional point-group token for
legacy workflows, but the authoritative GICForge report is now the group
computed by `symm.f`.

## Initial GIC Selection

The initial Fortran GIC generation follows the Merlino4 rule used by the
Python implementation. Bond stretches are not combined and are not reduced:
each bonded pair is kept as a primitive `R(i,j)` coordinate. Redundant
candidate sets are built for valence bends, linear bends, torsions,
out-of-plane terms and ring-specific coordinates, then reduced before output.
The final rank check is done against the actual B rows. This mirrors the
Python rule that comparable non-stretch coordinates must be grouped before
pruning.

The Python implementation still has richer signatures for grouping
coordinates: atom classes, ring tags, hydrogen tags, angle bins and local
topology are all available there. GICForge deliberately uses a conservative
Fortran77-compatible subset based on coordinate type plus atomic-number
patterns. This avoids false mixing while leaving the Python layer responsible
for higher-level topology policy.

## Type-Local Symmetrization

With `SYMMALL`, GICForge calls `SymOneGICBlock` before `OrdRed` only for
non-stretch coordinate families:

- bend
- linear bend
- torsion

The implementation is in `fortran/gicforge/gic_type_symmetry.f` and is called
from `fortran/gicforge/dina25.f`. It only acts on generic one-term coordinates
and only groups coordinates with the same family and compatible atomic-number
signature. Specialized ring/puckering/butterfly coordinates are left unchanged.
Out-of-plane signatures are implemented in the shared symmetry routine, but the
driver does not invoke them yet because `PrtOut` still cannot serialize OOP
linear combinations to Gaussian input.

Bond stretches are intentionally excluded from this local SALC step. They are
written as primitive one-term `R(i,j)` coordinates in the raw `gauin`; the later
deterministic Python post-check assigns irreducible representations and, when
required for Gaussian optimization input, writes the symmetry-adapted
`gauin.symm`.

For each homogeneous non-stretch group, the first GIC becomes the normalized
symmetric sum and the following GICs become orthonormal adjacent differences.
`OrdRed` and the later B-rank pruning then remove any remaining dependencies.
This preserves the required rule: no symmetrization or pruning step combines
different coordinate types.

## Residual Redundancy Pruning

After the ordinary GIC construction, GICForge calls `PruneGICBlocks` before
writing Gaussian input. The routine builds the B matrix for the current GIC
candidates and applies a modified Gram-Schmidt rank test separately to each
non-stretch coordinate family:

- bend
- linear bend
- torsion, including ring puckering and butterfly coordinates
- out-of-plane

Stretch rows are reported as kept primitive coordinates and are not pruned.
Dependent non-stretch rows are removed from their own family and the remaining
arrays are compacted in place. This is intentionally analogous to the Python
`prune_mode=svd/g` policy, but preserves strict Fortran77 implementation and
does not combine heterogeneous coordinate types.

## Python Symmetry Contract

After the Fortran77 backend has written `gauin`, Merlino runs a deterministic
GICForge post-check that assigns irreducible representations and writes:

- `gauin.symm`
- `gicsym`
- `gic_symmetry_diagnostics.json`

The number of final coordinates in each irrep is not chosen by the
semiexperimental solver. It is computed first from the vibrational
representation:

```text
Gamma_vib = Gamma_3N - Gamma_trans - Gamma_rot
```

For a non-linear molecule the sum of those counts must be `3N-6`. The
symmetry post-check also preserves the final GICForge coordinate-class counts:
the number of stretches, bends, linear bends, torsions and out-of-plane
coordinates in `gauin.symm` must match the non-symmetrized GICForge basis.
This prevents one class, for example stretches, from filling the rank that
belongs to torsions. If either the irrep counts or the coordinate-class counts
cannot be reached, the run fails. For Gaussian optimization, `gauin.symm`
always writes the totally symmetric coordinates first, then a blank separator,
then all other coordinates; the semiexperimental fit and geometry optimization
use only the GICForge-assigned totally symmetric block.

The post-check is intentionally reproducible:

- symmetry operations are sorted in a canonical order before projection
- numerical tolerances are centralized in `merlino_gic/gic_symmetry.py`
- repeated runs on the same `gauin`/`xyzin` write byte-identical `gauin.symm`,
  `gicsym` and diagnostics
- `gicsym` records the deterministic source block used for each coordinate

The projection hierarchy is deliberately restricted. A coordinate is written
only if it can be represented by one of these blocks:

- direct primitive projection from the GICForge expression
- same coordinate type (`bond`, `angle`, `linear_bend`, `dihedral`,
  `out_of_plane`)
- for non-totally-symmetric species only, an endocyclic ring block containing
  ring valence angles, ring dihedrals and ring-related improper/OOP
  coordinates
- for non-totally-symmetric species only, a local out-of-plane block containing
  dihedral-improper and OOP coordinates around the same center

There is no global least-squares fallback over all primitives. If these blocks
cannot generate the theoretical irrep counts, GICForge fails and the backend
must be fixed before the workflow proceeds.

`gic_symmetry_diagnostics.json` contains `strict_clean`, the irrep targets,
the obtained irrep counts, the coordinate-class targets, the obtained
coordinate-class counts and the source-block counts. `strict_clean=true` means
no unsupported/global reconstruction was used; it does not require every
coordinate to be a one-step primitive permutation, because type-local
Cartesian projection is the deterministic way to validate the B-row symmetry.
