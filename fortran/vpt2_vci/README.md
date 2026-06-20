Fortran77 VPT2/VCI Core
=======================

This directory contains new standalone Fortran77 kernels for the Merlino4
VPT2/VCI backend. They are not copied from GDV.

- `gf_core.f`: small symmetric Jacobi kernel and a Wilson-GF helper for already
  independent coordinates.
- `vci_core.f`: product-basis generation and dense small-space VCI helpers.

Gaussian parsing, coordinate construction and tensor normalization are handled
by Python. The Fortran code receives only numerical arrays and will later share
the same contracts as the Python implementation.

Use:

```bash
cd fortran/vpt2_vci
./compile_check
```
