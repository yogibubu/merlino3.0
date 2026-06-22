Merlino 2.1 (Freeze documentation)
# Merlino 2.1 — Design notes and freeze documentation

This document summarizes the main architectural and conceptual changes
introduced in Merlino 2.1 up to the current freeze point.

It is meant as a temporary design document and will be integrated into
the main documentation at a later stage.

1. General philosophy
## General philosophy

Merlino 2.1 is designed around a strict separation between:
- scientific pipelines
- project orchestration
- graphical user interface (GUI)

Scientific code is developed and validated independently from the GUI.
The GUI acts only as a trigger and visualization layer.

2. The xyzin enriched molecular container
## The xyzin enriched molecular container

`xyzin` is the central molecular container in Merlino 2.1.

It follows a sectioned enrichment model:
- input geometry and #BASIC are preserved
- enrichment sections are regenerated on each run
- previous enrichment blocks are replaced to keep the file consistent

`xyzin` can be generated from different input types:
- SMILES
- XYZ
- Z-matrix
- quantum chemistry outputs

Typical sections
Typical sections include:

- XYZ geometry block
- `#BASIC` (charge, multiplicity, symmetry)
- `#SMILES` (input SMILES, if provided)
- `#ROTATIONAL` (rotational constants and rotor classification)
- `#THERMO` (translational + rotational + vibrational + total thermodynamics)
- `#TOPOLOGY` (topological analysis)
- `#ISOTOPOLOGUES` (parent/isotopologue definitions, optionally enriched with
  rotational constants, DeltaVib/DeltaEl corrections and experimental sigmas)
- `#SMILES_Synthons` (SMILES reconstructed from topology)

The presence of a section depends on the input type and the executed pipelines.

3. Pipeline-based architecture
## Pipeline-based architecture

Each scientific domain in Merlino 2.1 is implemented as an independent pipeline.

A pipeline:
- reads the current `xyzin`
- performs a specific analysis
- appends new sections to `xyzin`
- optionally writes a human-readable `.report` file

Pipelines are:
- fully testable outside the GUI
- reusable in CLI, batch, or GUI mode
- independent from each other

Current pipelines
Currently implemented and validated pipelines:

- Geometry pipeline
  - rotational constants
  - rotor classification
  - thermodynamics (translational, rotational, vibrational, total)
  - output: `#ROTATIONAL`, `#VIBRATIONAL`, `#THERMO`, `rotational.report`, `thermo.report`

- Topology pipeline
  - atomic descriptors
  - bonding pattern
  - ring detection
  - aromaticity
  - topology-derived SMILES
  - output: `#TOPOLOGY`, `#SMILES_Synthons`, `topology.report`

4. ProjectManager
## ProjectManager

The ProjectManager orchestrates the execution of scientific pipelines.

Responsibilities:
- react to a validated `xyzin`
- trigger the appropriate pipelines
- manage paths and working directories

The ProjectManager:
- does not perform calculations
- does not depend on GUI libraries
- does not parse scientific results

5. GUI responsibilities
## GUI responsibilities

The GUI:
- collects user input
- validates and updates `xyzin`
- triggers the ProjectManager
- provides minimal user feedback

The GUI does not:
- perform scientific calculations
- duplicate pipeline logic
- parse numerical results

After pipeline execution, the GUI displays a minimal confirmation dialog:
- "Rotational analysis: OK"
- "Thermo analysis: OK"
- "Topology analysis: OK"

6. Reports
## Reports

Some pipelines generate additional human-readable reports:

- `rotational.report`
- `thermo.report`
- `topology.report`

These files are written in the working directory and can be displayed
directly in the GUI without additional parsing.

7. Gaussian integration (status)
## Gaussian integration (status)

Gaussian integration is intentionally postponed.

Reasons:
- Gaussian requires additional user-defined keywords
- GUI extensions are necessary
- design decisions are needed before implementation

The current architecture is designed to support Gaussian integration
as an additional pipeline without modifying existing pipelines.

8. Freeze status
## Freeze status

At this stage, the following components are considered stable:

- `xyzin` format and semantics
- geometry pipeline
- topology pipeline
- ProjectManager orchestration
- GUI-triggered execution and feedback

## Merlino modules

Merlino is organized into three high-level modules:

- **Merlino/Topology**  
  Handles molecular representation and automatic structural analysis.
  This module includes geometry, topology, and `xyzin` enrichment and
  represents the current stable core of Merlino 2.1.

- **Merlino/QM**  
  Handles preparation and execution of quantum-mechanical calculations.
  This module will include advanced preparation steps and support multiple
  quantum chemistry backends.

- **Merlino/Spectra**  
  Handles spectral simulation and post-processing starting from QM outputs.
  This module will focus on the construction of observable spectra and
  comparison with experiments.
