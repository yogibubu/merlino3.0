Architecture — Merlino 3.0
Overview

Merlino 3.0 is designed around a strictly modular, file-based architecture whose primary goal is:

robust, reproducible preparation of molecular input data

The architecture deliberately avoids:

shared mutable state

hidden dependencies

implicit data flow

All components interact only through files, with a single authoritative data channel.

Core architectural principles
1. Single source of truth

The only shared data structure in Merlino 3.0 is the file:

working/xyzin


All modules:

read from xyzin

write to xyzin

never exchange molecular data directly in memory

This guarantees:

transparency

debuggability

reproducibility

2. Fixed working directory

The working directory is always:

<project root>/working


It is:

created at startup if missing

cleaned at the end of a session unless DEBUG is enabled

never user-selectable

This removes ambiguity and simplifies all components.

3. Directory-level modularity

Each top-level directory is self-contained and can be developed and tested independently.

merlino3.0/
├── app.py
├── manager.py
├── working/
├── gui/
├── geometry/
├── topology/
├── bin/
└── doc/


There are no sub-packages with hidden cross-dependencies.

Module responsibilities
app.py

entry point

activates the environment

delegates control to manager.py

manager.py

initializes the working directory

creates the initial xyzin

launches the GUI

performs final cleanup (unless DEBUG)

manager.py owns the lifecycle, not the data.

gui/

Responsible for all user interaction.

Sub-components:

MainWindow: confirmation and rollback logic

InputPanel: input selection, preview, commit

ViewerPanel: visualization

The GUI:

never parses chemistry files directly

never stores molecular data internally

never bypasses xyzin

geometry/

Reserved for:

geometry manipulation

coordinate transformations

structure analysis

Currently independent from GUI and readers.

topology/

Reserved for:

connectivity analysis

bond perception

graph-based molecular operations

Interacts only via files or explicit calls.

bin/

Holds:

external executables

wrappers

helper scripts

Not imported directly by Python modules.

doc/

Project documentation:

Architecture.md

GUI.md

Testing.md

Documentation reflects current reality only.

Readers architecture

Readers are responsible for converting external representations into xyzin.

Characteristics:

one reader per input type

stateless

no GUI awareness

no confirmation logic

Readers:

overwrite working/xyzin on success

raise exceptions on failure

never ask user input

GUI and readers interaction
User action
   ↓
InputPanel
   ↓
Reader
   ↓
working/xyzin (written)
   ↓
MainWindow confirmation
   ↓
Accept → Viewer refresh
Reject → rollback


This pipeline is strictly enforced.

Preview vs commit

Preview operations:

are non-destructive

never write to xyzin

never trigger confirmation dialogs

Commit operations:

always write to xyzin

always require confirmation

always support rollback

This distinction is fundamental and must not be violated.

Testing philosophy

automated tests cover:

readers

xyzin correctness

GUI logic is tested manually

no attempt is made to simulate Qt events automatically

This ensures:

fast test suite

stable CI

focus on scientific correctness

What is explicitly out of scope

Merlino 3.0 does not include:

project management

job submission

workflow automation

database storage

background calculations

These are intentionally excluded to keep the architecture minimal and robust.

Final remarks

Merlino 3.0 favors:

clarity over flexibility

explicitness over convenience

files over shared state

stability over features

The architecture is intentionally conservative to support long-term scientific work.
