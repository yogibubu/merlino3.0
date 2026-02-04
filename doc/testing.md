Testing Merlino 3.0
Overview

Merlino 3.0 uses pytest to validate the core functionality of the codebase, with particular emphasis on:

correctness of readers

correct generation and update of the xyzin file

consistency of charge and multiplicity handling

robustness against different input formats

At this stage, testing focuses on non-GUI logic, with a small set of
offscreen GUI unit tests.

Project assumptions (important)

The working directory is always:

<project root>/working


The only communication channel between modules is the file:

working/xyzin


Readers always overwrite working/xyzin on successful commit.

GUI logic (confirmation, rollback) is not part of automated tests.

Test layout
merlino3.0/
├── gui/
│   ├── test_readers.py
│   ├── test_gaussian.py
│   └── tests/
│       └── gaussian/
│           ├── h2o.out
│           ├── h2o.log
│           ├── h2o.fchk
│           └── ...
└── doc/
    └── Testing.md

Test files

gui/test_readers.py

Tests basic readers (e.g. XYZ)

Verifies correct overwrite of working/xyzin

Uses tmp_path to simulate project root + working directory

gui/test_gaussian.py

Tests Gaussian reader

Automatically scans all files in gui/tests/gaussian

Supported formats:

.out

.log

.fchk

.fch

Verifies:

XYZ geometry is written

#BASIC section is present

charge and multiplicity are correctly extracted

How to run tests

All tests are run from the project root:

pytest gui/test_readers.py
pytest gui/test_gaussian.py


Or all GUI-related tests:

pytest gui

Offscreen GUI tests (pytest-qt)

QT_QPA_PLATFORM=offscreen pytest gui/tests

Notes on Gaussian tests

If multiple files with the same stem exist (e.g. h2o.out and h2o.fchk):

each file is tested independently

the reader decides internally how to parse the file

the test only checks the final result written to xyzin

The test does not enforce priority between formats; it verifies correctness.

GUI testing

The GUI is tested both via pytest-qt (offscreen) and manually by running:

python app.py


and verifying:

correct SMILES preview

commit only on ENTER

correct confirmation dialog

rollback on rejection

correct handling of file-based inputs

This is a deliberate design choice to keep tests:

deterministic

fast

independent of Qt event handling

Adding new reader tests

To add a new reader test (e.g. Molpro, MRCC):

Create a directory under:

gui/tests/<program>/


Add representative output files

Copy and adapt test_gaussian.py

Ensure:

geometry is written

#BASIC is updated consistently

Final note

Testing in Merlino 3.0 is intentionally minimal, focused, and realistic.

The goal is not exhaustive coverage, but:

early detection of structural regressions

guaranteed correctness of xyzin

confidence in reader integration
