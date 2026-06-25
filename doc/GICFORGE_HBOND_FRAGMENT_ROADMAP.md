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

## Added Reference: Skeletal Natural Internal Coordinates

The local reference
`doc/papers/references/berces_1998_metal_complex_natural_internal_coordinates.pdf`
is now part of the workspace references.  Its relevant lesson for GICForge is
not limited to transition-metal complexes.  Berces shows that the problematic
coordinates are the skeletal degrees of freedom: the translation and rotation
of a ligand or subfragment relative to another molecular object.  These should
be represented by internal coordinates built from a geometrical reference point
of the bonded atoms, not by a dense set of ordinary atom-atom pseudo-bonds.

The design consequences for GICForge are:

- a fragment or subfragment may be an internal object of a connected molecule,
  not only a disconnected molecule;
- non-atomic centers are part of the coordinate definition, but they must not
  become dummy atoms in the Merlino topology;
- geometrical centers are preferred over centers of mass, so isotopic
  substitution does not change the coordinate definitions;
- the center atoms are selected by the interaction mode: for an eta-n ligand,
  use the n atoms actually bound or coordinated, not necessarily all atoms in
  the ligand;
- skeletal motion is a structured six-dimensional block whenever both partners
  are nonlinear: three translational-type coordinates and three
  rotational-type coordinates;
- centroid/frame coordinates should be preferred to sums of many metal-ligand
  or contact pseudo-bonds, because the latter couple strongly to intra-fragment
  stretches and bends;
- local-symmetry combinations of equivalent angles or torsions are a useful
  starting point, but the final implementation must still pass the GICForge
  rank and symmetry checks.

## Decision After geomeTRIC/TRIC Review

geomeTRIC's TRIC model is the right reference point for disconnected systems:
it treats fragments as collective objects and adds translations plus rotations
for each non-reference fragment.  In the TRIC paper and documentation, the
translation coordinates are fragment centroids and the rotations are represented
with quaternion exponential-map coordinates.  This is robust for molecular
clusters and non-covalent optimizations, but it is not the most chemically
specific representation for every Merlino use case.

GICForge should keep the useful TRIC idea, but generalize it in a more local
and chemically explicit way:

- a "fragment" may be a whole disconnected molecule, a covalent subfragment, a
  bond-centered object, a ring-centered object, or another user/topology-defined
  local object;
- each fragment/subfragment has a deterministic center and an orientation frame;
- the frame is defined by the center plus two orientation atoms whenever
  possible, so centers at a bond midpoint, ring centroid, coordination center,
  lone-pair direction or other non-atomic site can be used without introducing
  dummy atoms into the molecular graph;
- inter-fragment coordinates are built between these local frames rather than
  only between disconnected-molecule centroids.

This differs from a pure TRIC policy.  TRIC is mainly a rigid-fragment
translation/rotation layer.  The GICForge extension should be a block-aware
local-frame coordinate system that can coexist with ordinary non-redundant GICs,
symmetry labels and Gaussian-readable coordinate definitions.

Berces' skeletal-coordinate construction is the complementary reference for
connected systems, haptic ligands, bridges and coordination centers.  In this
case the "inter-fragment" coordinate is often internal to a single covalent
component.  The correct object is therefore not a disconnected-fragment TRIC
coordinate, but a skeletal local-frame coordinate between a center atom or
frame and a subfragment centroid/frame.

## Guiding Principle

H-bonds and other weak contacts must not silently become ordinary covalent
bonds. They should be stored as non-covalent contact annotations and used
explicitly by modules that request them.  When requested for GICForge, they
enter as lower-weight contact coordinates and pseudo-cycle descriptors in their
own block, not as ordinary covalent bonds.

The coordinate model should have separate blocks:

- `intra_fragment`: covalent GICs within each connected covalent fragment.
- `intra_subfragment_skeletal`: skeletal coordinates inside one covalent
  component, for example metal-ligand eta-n motion, bridge motion relative to a
  ring, substituent motion relative to a bond center, or ligand motion relative
  to a coordination center.
- `local_frame`: definitions of user/topology fragments and subfragments, each
  with center, two orientation atoms and fallback rules for degenerate frames.
- `inter_frame`: relative distances, angular coordinates and torsions between
  local frames or between atoms and local frames.
- `inter_fragment_rigid`: TRIC-like relative translations and rotations for
  truly disconnected fragments.
- `inter_fragment_contacts`: directed pseudo-bonds and angular descriptors for
  H-bonds and other chemically meaningful contacts.
- `pseudo_cycle_contacts`: H-bond/weak-contact pseudo-ring closure coordinates,
  kept below covalent GICs in priority and metric weight.
- `constraints_or_predicates`: optional contact-derived relations that guide a
  fit or geometry correction without entering the active GIC space as ordinary
  covalent coordinates.

Redundancy removal and symmetry adaptation must be block-aware. Stretching
coordinates must not replace torsions, and contact coordinates must not replace
intra-fragment coordinates.

## Local Frame Model

Recommended frame record:

```text
type = fragment_frame
id = F1
atoms = 1,2,3,4,5,6
center = atom:3 | bond:3-4 | ring:1,2,3,4,5,6 | centroid:...
axis_atom_1 = 1
axis_atom_2 = 2
orientation_policy = center_axis_plane
weight = 1.0
source = user | topology | ring | bond | coordination
hapticity = none | eta1 | eta2 | eta3 | eta4 | eta5 | eta6 | ...
center_atoms = 1,2,3
skeleton_partner = atom:7 | frame:F2 | none
```

Frame construction:

- origin `O`: atom position, bond midpoint, ring centroid, coordination-center
  centroid, or user-defined weighted centroid;
- default centroid policy: geometrical center of the atoms defining the
  interaction or coordination mode, not center of mass;
- first axis: normalized vector from `O` to `axis_atom_1`;
- second direction: component of `axis_atom_2 - O` orthogonal to the first
  axis;
- third axis: right-handed cross product;
- if the two orientation atoms are nearly collinear with the center, fall back
  deterministically to the best-conditioned atom pair in the fragment and
  report the substitution.

Initial coordinate families:

- frame--frame distance `R(F1,F2)`;
- frame--frame orientation angles, preferably expressed as stable
  quaternion/exponential-map increments for large rotations;
- skeleton stretch between a center atom/frame and a subfragment centroid;
- skeleton bend/tilt generated from equivalent center-centroid-atom angles and
  then reduced into rank-valid local combinations;
- skeleton internal rotation from equivalent atom-centroid-center-reference
  torsions;
- skeleton out-of-plane motion when a centroid leaves a coordination plane;
- atom--frame distance/angle/dihedral for coordination to a bond center, ring
  center or metal/coordination center;
- optional frame-local Cartesian displacements for SYCART-like workflows.

The frame definition is metadata, not an extra atom in the covalent graph.
Dummy atoms may be written to Gaussian only as an output representation if that
is the cleanest way to express a GIC, but they are not part of Merlino topology.

## Automatic Detection of Non-Atomic Centers

User-defined frames must always override automatic perception.  When no
explicit `#FRAMES` or `#SKELETON` section is present, GICForge should generate
candidate bond, ring, coordination and stacking frames from the covalent graph
and the current Cartesian geometry.  Automatic frames are candidates only: they
enter the active coordinate model only if requested by keyword and if their
block passes the rank, symmetry and priority checks.

Recommended candidate record:

```text
type = auto_frame_candidate
kind = bond_center | ring_center | coordination_center | eta_center | stacking
center_atoms = ...
partner = atom:i | bond:i-j | ring:... | frame:F
score = 0.0 .. 1.0
distance_spread
plane_rms
normal_alignment
source = geometry+topology
```

### Bond-Center Candidates

A bond-centered frame should be proposed when an external atom or frame sees
the two atoms of a bond almost equivalently.

Detection:

- start from topological bonds `i-j`;
- for each external atom or frame center `X`, compute `d(X,i)` and `d(X,j)`;
- accept a bond-center candidate if `abs(d(X,i)-d(X,j))` is below a relative
  tolerance, for example 5-8% of the mean distance, and `X` is not simply a
  normal covalent neighbor of only one endpoint;
- prefer bonds whose midpoint lies close to the line or approach direction from
  `X`;
- reject candidates where the midpoint frame is nearly collinear with all
  possible orientation atoms.

Typical uses:

- eta2 coordination;
- substituent or bridge orientation about a bond center;
- approach of a weak contact to a multiple bond;
- internal rotors where the chemically meaningful center is a bond midpoint.

### Ring-Center and Eta-n Candidates

A ring-centered frame should be proposed when an atom or frame is nearly
equidistant from several atoms belonging to the same perceived ring or fused
ring face.  The center atoms must be the atoms participating in the interaction,
not necessarily every atom in the ring.

Detection:

- use the existing ring perception to list simple rings and fused-ring faces;
- for each external atom or frame `X`, compute distances to ring atoms;
- search subsets of size 3-8, with priority to chemically common eta2, eta3,
  eta4, eta5 and eta6 subsets;
- accept a subset if the distance spread is small, for example max-min below
  8-12% of the mean distance, and the subset is geometrically coherent;
- require the projection of `X` onto the ring/subset plane to fall inside or
  close to the polygon/convex hull unless the user explicitly allows edge
  approaches;
- for planar rings, use the ring normal and in-plane principal axes as the
  default frame; for non-planar rings, use the least-squares plane only if the
  plane RMS is below a threshold, otherwise split into smaller local subsets;
- for fused PAHs, evaluate each smallest cycle and each fused face, then keep
  the highest-scoring non-overlapping candidates.

The scoring should combine distance equivalence, plane quality, projection
inside the face, and chemical plausibility of the partner atom.  This lets
GICForge detect eta-n coordination and ring-centered interactions without
hardcoding specific ligands.

### Coordination-Center Candidates

Coordination-centered frames are needed for metal atoms, hypercoordinate
centers, centers above coordination number four, and other cases where ordinary
valence-angle templates are insufficient.

Detection:

- start from atoms with high coordination number, metal-like element classes,
  or user-requested `SKELETALGIC`;
- group neighbors by distance shell and element/type compatibility;
- when a neighbor subset has nearly equal center-neighbor distances, propose a
  local coordination frame based on that subset;
- for eta-n ligands, treat the ligand centroid/frame as one coordination
  object rather than expanding it into many independent pseudo-bonds;
- keep monohapto eta1 cases as ordinary atom-atom coordination unless a
  frame-based representation gives a better block rank.

### Automatic Stacking Frames

Stacking should be detected as an inter-frame interaction between two ring or
pi-system frames, not as a collection of atom-atom contacts.

Detection:

- build ring/pi frames for aromatic or conjugated cycles and for user-defined
  pi fragments;
- compare every pair of frames in different fragments or distant subfragments;
- accept a stacking candidate when the centroid distance is in a plausible
  non-covalent range, the normals are near parallel or antiparallel, and the
  projected centroid offset is compatible with face-to-face or slipped stacking;
- classify the geometry as face-to-face, slipped-parallel or edge-to-face using
  normal alignment and centroid-offset components;
- generate frame-frame distance, lateral offset and relative rotation
  coordinates as a low-priority inter-frame block;
- do not create covalent rings, covalent torsions or H-bond pseudo-cycles from
  stacking unless a separate contact rule requests them.

For PAH stacks and benzene-dimer-like systems, the preferred primitives are
therefore ring-center distance, two lateral displacement coordinates and one
relative twist/tilt block, with redundancy pruning performed inside the
stacking block.

### Priority and False-Positive Control

Automatic frame generation must be conservative:

- explicit user frames outrank automatic frames;
- covalent intra-fragment GICs outrank skeletal, stacking and contact frames;
- connected skeletal frames outrank weak-contact pseudo-cycles;
- weak-contact and stacking frames are inactive unless requested by keyword;
- candidates below a score threshold are reported in the manifest but not used;
- overlapping candidates are resolved by score, smaller chemically specific
  subsets before large generic subsets, and deterministic atom-order tie breaks;
- every accepted candidate must report center atoms, partner, score and reason
  in `provout`, `gauin` comments when relevant, and the JSON manifest.

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

These contacts can be used for BDPCS3 geometry correction immediately.  For
GICForge they should enter only through an explicit keyword, for example
`CONTACTGIC` or `HBONDCYC`, and then only in the contact/pseudo-cycle block.

The intended H-bond coordinate policy is:

- H...A stretch: included as a lower-weight pseudo-bond, not as a covalent
  stretch;
- D-H...A valence angle: included as a directed donor-angle descriptor;
- H...A-X or D-H...A-X torsion/dihedral: included when the donor and acceptor
  neighborhoods define a stable pseudo-ring closure;
- contact metric weight: lower than covalent stretches, valence bends and
  torsions; smooth distance and angle damping from the common parameter table;
- pseudo-cycle closure: allowed only in `pseudo_cycle_contacts`, never in the
  covalent ring finder.

## Fragment and Frame Coordinates

For disconnected systems, the robust default should remain TRIC-like:

- choose a reference fragment deterministically, usually the largest fragment
  with canonical atom-order tie breaks;
- add three relative translations for each non-reference fragment;
- add three relative rotations using the quaternion exponential map;
- keep these coordinates in an `inter_fragment_rigid` block.

For covalently connected systems, frame coordinates should be more general than
TRIC:

- bond-centered frames for substituent orientation, ligand approach, bridge
  coordinates and internal rotors;
- ring-centered frames for pi-stacking, ring puckering couplings and
  ring-substituent orientation;
- coordination-centered frames for metal centers, hypervalent centers and
  coordination numbers above four;
- hapticity-centered frames for eta-n ligands, where the center is defined by
  the atoms actually bound to the metal or coordination site;
- user-defined subfragment frames for constrained fitting or geometry
  optimization.

The default connected-fragment policy should follow the skeletal NIC model:

- build ordinary intra-fragment GICs inside each ligand/subfragment first;
- build skeletal stretch, tilt, internal rotation and skeletal bend/out-of-plane
  coordinates between the subfragment frame and the partner center/frame;
- keep skeletal coordinates in their own block, after covalent intra-fragment
  coordinates and before weak-contact pseudo-cycles;
- use local combinations of symmetry-equivalent primitive angles/torsions when
  the topology gives equivalent atoms, but let the existing rank and symmetry
  machinery decide the final independent rows;
- avoid replacing a centroid/frame skeletal coordinate by many atom-center
  pseudo-bonds except as an explicitly requested diagnostic.

Pseudo-bond contacts and frame/quaternion coordinates solve different problems:

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

#FRAMES
SCHEMA merlino.xyzin.frames.v1
FRAME ring1 ATOMS=1,2,3,4,5,6 CENTER=ring:1,2,3,4,5,6 ORIENT=1,2
FRAME bond1 ATOMS=3,4,7,8 CENTER=bond:3-4 ORIENT=7,8
FRAME coord1 ATOMS=1,2,3,4,5 CENTER=atom:1 ORIENT=2,3
FRAME eta5_cp ATOMS=2,3,4,5,6,7,8,9,10,11 CENTER=centroid:2,3,4,5,6 ORIENT=2,3 HAPTICITY=eta5 PARTNER=atom:1

#CONTACTS
SCHEMA merlino.xyzin.contacts.v1
HBOND donor=2 hydrogen=3 acceptor=5 weight=0.83 distance_HA=1.92 angle_DHA=168.4

#SKELETON
SCHEMA merlino.xyzin.skeleton.v1
SKELETON center=atom:1 frame=eta5_cp mode=eta5 COORDS=stretch,tilt,rotation
```

GICForge should read these sections when present, or generate them from the
same shared topology/contact library when absent.

## TODO / Implementation Phases

1. **Policy and schema**
   - Rename this roadmap into a formal GICForge design note once the policy is
     stable.
   - Define `FragmentRecord`, `FrameRecord`, `SkeletonRecord`,
     `ContactRecord` and `PseudoCycleRecord` schemas.
   - Add optional `#FRAMES` to `xyzin` without breaking old inputs.
   - Add optional `#SKELETON` to `xyzin` for connected subfragment skeletal
     coordinates.
   - Decide keyword names:
     `FRAMES` for local frame generation, `CONTACTGIC` for contact coordinates,
     `SKELETALGIC` for centroid/frame skeletal coordinates, `HBONDCYC` for
     H-bond pseudo-cycle closure, and `NOFRAGTRIC` for disabling
     disconnected-fragment rigid coordinates.

2. **Common topology/contact library**
   - Move fragment detection, H-bond detection, and smooth contact weights into
     `merlino_core`.
   - Make Python GICForge, BDPCS3, BSR, and GUI modules call this one library.
   - Export the same parameters to Fortran77 include files.
   - Add frame builders for atom, bond, ring, centroid and user-defined centers.
   - Add hapticity-aware center builders that use the coordinated atoms rather
     than the full ligand by default.
   - Add automatic candidate detectors for bond centers, ring/eta-n centers,
     coordination centers and stacking frames, with scores and rejection
     reasons.
   - Add deterministic orientation-atom selection and collinearity fallback.
   - Add a chain-rule derivative path for non-atomic centers so B rows are
     distributed to real atoms only.

3. **Python GICForge reference behavior**
   - Add block labels to every primitive and final GIC:
     `scope`, `coordinate_type`, `fragment_ids`, `contact_id`.
   - Keep intra-fragment, intra-subfragment skeletal, inter-fragment rigid, and
     contact blocks separate through redundancy removal and symmetry adaptation.
   - Add frame coordinates as first-class primitive kinds:
     `frame_distance`, `frame_angle`, `frame_torsion`, `frame_rotation`,
     `atom_frame_distance`, `atom_frame_angle`, `atom_frame_torsion`.
   - Add skeletal primitive kinds:
     `skeleton_stretch`, `skeleton_tilt`, `skeleton_rotation`,
     `skeleton_bend`, `skeleton_out_of_plane`.
   - Prefer centroid/frame skeletal rows over atom-center pseudo-bond sums; keep
     pseudo-bond representations only as diagnostics or explicit fallback.
   - Add lower metric weights for H-bond/contact primitives.
   - Implement pseudo-cycle contact generation for D-H...A plus donor/acceptor
     angular and torsional closure descriptors.
   - Add deterministic tests for fragment count, contact count, block ranks,
     total rank, and symmetry labels.
   - Add tests that distinguish true eta/ring-center candidates from ordinary
     atom-atom contacts and from false positives caused by accidental similar
     distances.

4. **Fortran77 parity**
   - Mirror the Python contact and fragment rules in strict Fortran77.
   - Keep Fortran data structures fixed-size and deterministic.
   - Add Python/Fortran comparison tests for primitive lists, block ranks, final
     GIC labels, and B-matrix rows.
   - Implement frame-center and orientation fallback with identical numerical
     tolerances.
   - Implement non-atomic center derivatives by chain rule and distribute every
     center derivative to the defining real atoms.
   - Keep geometrical-centroid definitions isotope-independent in both Python
     and Fortran.
   - Ensure `provout`, `gauin`, `bmat.out` and the manifest carry frame/contact
     block labels.

5. **Gaussian and downstream interfaces**
   - Write Gaussian GIC input with clear inactive/active partitioning:
     intra-fragment active coordinates first, then optional inter-fragment
     coordinates, then inactive diagnostics/contact descriptors.
   - Preserve symmetry labels and block labels in readable reports and JSON.
   - Decide whether non-atomic centers are written as Gaussian dummy atoms or as
     explicit GIC expressions, keeping Merlino topology dummy-free.
   - Make BSR and future geometry optimizers request frame/contact blocks
     explicitly rather than receiving them by default.

6. **Validation set**
   - Water dimer: pure H-bonded two-fragment case.
   - Formic acid--water or ammonia--water: strong directed H-bond.
   - Glycine conformers: intramolecular H-bond/contact without fragment merge.
   - A van der Waals dimer: inter-fragment rigid coordinates without H-bond.
   - A bridged/multifragment edge case to test rank stability.
   - Benzene dimer or stacked PAH pair: ring-center frames and pi-stacking
     orientation.
   - Parallel, slipped and edge-to-face benzene dimers: stacking classification
     and low-priority inter-frame block ranks.
   - Metal/ligand or hypercoordinate model: coordination-center frame.
   - Ferrocene or benzene-chromium-tricarbonyl style eta-n ligand: skeletal
     stretch, tilt and internal rotation against centroid/frame alternatives.
   - Biphenyl or norbornane derivative: bond-centered frame and bridge
     orientation.
   - A deliberately near-linear H-bond case to test angle/torsion stability.

## Non-Negotiable Rules

- Python and Fortran77 GICForge must produce identical coordinate definitions.
- H-bonds do not alter covalent fragment identity.  If the user explicitly
  requests pseudo-cycle behavior, it is generated only in the contact block.
- Redundancy removal is block-aware and cannot substitute one coordinate type
  for another.
- Symmetry labels must be assigned to all final coordinates, not only A1.
- BSR and Gaussian optimizations may select only the totally symmetric active
  subset, but reports must show all coordinates and their symmetries.
- Frame centers and orientation atoms are deterministic and appear in the
  manifest.
- Automatically detected centers are candidates with scores and rejection
  reasons; they are not silently added to the active space.
- Non-atomic centers are geometrical centers unless the user explicitly requests
  mass weighting; default GIC definitions must be isotope-independent.
- Eta-n and coordination-center fragments use the atoms actually involved in
  the interaction to define the centroid/frame.
- Skeletal local-frame coordinates are preferred over dense pseudo-bond
  representations for connected fragment motion.
- Non-atomic centers are never topological dummy atoms.  Their B-matrix
  derivatives are distributed to real atoms by the chain rule.
- Contact coordinates have lower priority and lower metric weight than
  covalent GICs unless the user explicitly overrides the weights.
