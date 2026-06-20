# GDV VPT2/VCI Source Map

Merlino4 uses the local GDV source checkout as the reference implementation for
the harmonic, anharmonic and vibrational-CI refactor.

Local source root:

```text
/Users/vincenzobarone/gdv_j32p/gdv
```

## Harmonic Internal Coordinates

- Imported active source: `fortran/harmonic_internal/gf.f`
- Original GDV source: `dinautil.F`
- Relevant decks:
  - `DNICGF`: builds internal-coordinate F/G matrices from Cartesian gradient,
    Cartesian Hessian, B matrix and B-matrix derivatives.
  - `DNICFq`: Wilson GF harmonic frequencies and internal-coordinate normal
    modes.

## Anharmonic Data And VPT2 Support

- Original GDV source: `dinautil.F`
- Relevant decks:
  - `AnhFIO`: anharmonic read/write file access for cubic/quartic force
    constants and Coriolis terms.
  - `VibAlpha`, `VibCNM`, `VibCor`, `VibINM`, `VPT2En`: vibro-rotational,
    Coriolis and VPT2 support routines.

## Vibrational CI

- Original GDV source: `l717.F`
- Relevant decks:
  - `VCI1MI`: monomode harmonic-oscillator integrals.
  - `VCIDrv`: VCI main driver.
  - `VCIGen`: automatic vibrational-state generation.
  - `VCIInt`: VCI Hamiltonian matrix elements.
  - `VCIPrt`: VCI output.
  - `VCIPT2`: mixed perturbative treatment.
  - `VCIVar`: pure variational treatment.
  - `VPTCrs`, `VPTRHS`, `VPTWFC`: VPT2 wavefunction coefficients.

Current limitation in `VCIDrv`: when the VCI space exceeds the dense threshold,
the code still stops with "Davidson not yet available". Merlino4 should replace
that branch with a standalone Davidson matrix-vector backend.

## Davidson/Subspace Diagonalization

- Original GDV source: `utilnz.F`
- Relevant deck:
  - `NHDiag`: non-Hermitian subspace diagonalization used by Davidson-like
    iterations elsewhere in GDV.

`utilnz.F` is about 26 MB and should not be copied wholesale into Merlino4.
Extract only the routines needed by the standalone Davidson backend and keep
their dependencies explicit.
