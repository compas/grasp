Radial grid in GRASP
====================

The grid is defined using the following parameters which are stored in the `grid_C` module in `libmod`.

* \\(R\\) (`RNT`): defines the scale of the grid. It has a dimension of length, and so its value is given in atomic units (Bohr radii). We'll denote it by \\(R\\) in the following equations.
* \\(H\\) (`H`):
* \\(H\\) (`H`):
* \\(H_{\mathrm{P}}\\) (`HP`): determines whether the grid is a simple exponential grid (\\(H_{\mathrm{P}} = 1\\)) or a more complex asymptotically linear-exponential grid (\\(H_{\mathrm{P}} > 0\\)).

* The grid parameters lives in `grid_C`
* The grid is also affected by compile-time parameters from  the `parameter_def` module.
  - `NNNP`:
  - `NNN1`: determines the sizes of the arrays, and should be set to `NNNP + 10`.

They live in the `grid_C` global module.

* The maximum number of

* The `RADGRD` subroutine from `lib9290`.

Arrays:

* `R`:
* `RP`:
* `RPOR`:

\\[
    r(i) = R \left[ \mathrm{e}^{(i-1)H} - 1 \right]
\\]

## Default grid parameters

## Tips and tricks

### A finer grid
