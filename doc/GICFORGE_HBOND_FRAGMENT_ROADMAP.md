# GICForge H-Bond and Fragment-Coordinate Roadmap

## Current Position

GICForge already has partial support for the next development step:

- the Python side already adds TRIC-style fragment coordinates for disconnected
  systems: fragment translations plus rotations from the quaternion exponential
  map;
- the Fortran77 GICForge code already detects fragments and has H-bond
  detection routines that can bridge fragments without changing the covalent
  topology by default;
- BDPCS3 H-bond parameters are centralized in `merlino_core.parameters` and
  exported to Fortran include files.

The next refactor must consolidate these pieces into a single molecular
coordinate model used by all Merlino modules.

## Guiding Principle

H-bonds and other weak contacts must not silently become ordinary covalent
bonds. They should be stored as non-covalent contact annotations and used
explicitly by modules that request them.

The coordinate model should have separate blocks:

- `intra_fragment`: covalent GICs within each connected covalent fragment.
- `inter_fragment_rigid`: relative translations and rotations between
  fragments, using TRIC/quaternion exp-map coordinates.
- `inter_fragment_contacts`: directed pseudo-bonds and angular descriptors for
  H-bonds and other chemically meaningful contacts.
- `constraints_or_predicates`: optional contact-derived relations that guide a
  fit or geometry correction without entering the active GIC space as ordinary
  covalent coordinates.

Redundancy removal and symmetry adaptation must be block-aware. Stretching
coordinates must not replace torsions, and contact coordinates must not replace
intra-fragment coordinates.

## H-Bond Model

H-bonds should be identified by the common topology/contact engine, not
independently in each module.

Recommended contact record:

```text
type = hbond
donor = D
hydrogen = H
acceptor = A
fragments = fD, fA
distance_HA
angle_DHA
weight_distance
weight_angle
weight_total
```

The weight should be smooth:

- distance switching from the shared BDPCS3 error-function model;
- angular damping with a threshold around 150 degrees for D-H...A;
- contact type correction from the centralized parameter table.

The present BDPCS3 policy remains:

- O-H...O uses the reoptimized correction value `-0.055`;
- N-H...N initially uses the same value;
- mixed N-H...O / O-H...N values use the average until reoptimized.

These contacts can be used for BDPCS3 geometry correction immediately, but
should only enter GICForge active coordinates through an explicit option.

## Fragment Coordinates

For disconnected systems, the robust default should be TRIC-like:

- choose a reference fragment deterministically, usually the largest fragment
  with canonical atom-order tie breaks;
- add three relative translations for each non-reference fragment;
- add three relative rotations using the quaternion exponential map;
- keep these coordinates in an `inter_fragment_rigid` block.

Pseudo-bond contacts and quaternion coordinates solve different problems:

- quaternions define the rigid relative placement of fragments and give the
  correct six degrees of freedom per free fragment;
- pseudo-bonds encode chemically directed interactions such as H-bonds,
  halogen bonds, ion pairs, or coordination contacts.

Therefore both should be available, but never mixed before rank reduction.

## Proposed `xyzin` Sections

The common Merlino container should be extended incrementally with optional
sections:

```text
#FRAGMENTS
SCHEMA merlino.xyzin.fragments.v1
FRAGMENT 1 ATOMS=1,2,3
FRAGMENT 2 ATOMS=4,5,6

#CONTACTS
SCHEMA merlino.xyzin.contacts.v1
HBOND donor=2 hydrogen=3 acceptor=5 weight=0.83 distance_HA=1.92 angle_DHA=168.4
```

GICForge should read these sections when present, or generate them from the
same shared topology/contact library when absent.

## Implementation Phases

1. **Repository cleanup and freeze**
   - Finish the current `xyzin`/isotopologue changes and commit them.
   - Keep ignored runtime data out of git.
   - Decide separately whether to untrack `working/xyzin`, because it is a
     runtime file but currently tracked.

2. **Common topology/contact library**
   - Move fragment detection, H-bond detection, and smooth contact weights into
     `merlino_core`.
   - Make Python GICForge, BDPCS3, BSR, and GUI modules call this one library.
   - Export the same parameters to Fortran77 include files.

3. **Python GICForge reference behavior**
   - Add block labels to every primitive and final GIC:
     `scope`, `coordinate_type`, `fragment_ids`, `contact_id`.
   - Keep intra-fragment, inter-fragment rigid, and contact blocks separate
     through redundancy removal and symmetry adaptation.
   - Add deterministic tests for fragment count, contact count, block ranks,
     total rank, and symmetry labels.

4. **Fortran77 parity**
   - Mirror the Python contact and fragment rules in strict Fortran77.
   - Keep Fortran data structures fixed-size and deterministic.
   - Add Python/Fortran comparison tests for primitive lists, block ranks, final
     GIC labels, and B-matrix rows.

5. **Gaussian and downstream interfaces**
   - Write Gaussian GIC input with clear inactive/active partitioning:
     intra-fragment active coordinates first, then optional inter-fragment
     coordinates, then inactive diagnostics/contact descriptors.
   - Preserve symmetry labels and block labels in readable reports and JSON.

6. **Validation set**
   - Water dimer: pure H-bonded two-fragment case.
   - Formic acid--water or ammonia--water: strong directed H-bond.
   - Glycine conformers: intramolecular H-bond/contact without fragment merge.
   - A van der Waals dimer: inter-fragment rigid coordinates without H-bond.
   - A bridged/multifragment edge case to test rank stability.

## Non-Negotiable Rules

- Python and Fortran77 GICForge must produce identical coordinate definitions.
- H-bonds do not alter covalent fragment identity unless the user explicitly
  requests pseudo-cycle behavior.
- Redundancy removal is block-aware and cannot substitute one coordinate type
  for another.
- Symmetry labels must be assigned to all final coordinates, not only A1.
- BSR and Gaussian optimizations may select only the totally symmetric active
  subset, but reports must show all coordinates and their symmetries.
