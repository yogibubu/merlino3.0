# BDPCS3 Fortran Backend

This backend is the Fortran77 mirror of the Python implementation in
`merlino_fit/survibfit/modify_geom.py`, specifically
`bdpcs3_delta_and_order_updated`.

The public routines are:

- `BDPCS3_UPDATED_BOND(ZI,ZJ,R,DELTA,BONDORD,RBD)`
- `BDPCS3_UPDATED_BOND_OVR(ZI,ZJ,R,BOOVR,IOVR,DELTA,BONDORD,RBD)`

`BDPCS3_UPDATED_BOND_OVR` reproduces the Python topology bond-order override:
when `IOVR` is nonzero, `BONDORD=max(BOOVR, exp((r_cov-r)/0.30))`.

Hydrogen-bond targets are included through:

- `BDPCS3_HBOND_TARGET(ZD,ZA,RHY,ANGXHY,DELTA,RBD)`

The target is the H...Y distance in X-H...Y. The current angular gate is
X-H-Y >= 150 deg and the correction is damped by an error function centered at
3.0 Angstrom so long contacts go continuously to zero. Current values are
O-H...O = -0.055 Angstrom and N-H...N = -0.055 Angstrom; mixed N/O cases use
their average.

The Fortran constants are in
`merlino_core/parameters/fortran/bdpcs3_hbond_params.inc`, which is kept
consistent with `merlino_core/parameters/bdpcs3_hbond.toml`. Regenerate the
include with:

```bash
python -m merlino_core.parameters.generate_fortran_includes
```

Build:

```bash
./compile_MAC
```
