# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    varioplot(γ; [options])

Plot the variogram or varioplane `γ` with given `options`.

## Empirical variogram options:

* `vcolor` - color of variogram
* `psize`  - size of points of variogram
* `ssize`  - size of segments of variogram
* `tshow`  - show text counts
* `tsize`  - size of text counts
* `hshow`  - show histogram
* `hcolor` - color of histogram

## Empirical varioplane options:

* `vscheme` - color scheme of varioplane
* `rshow`   - show range of theoretical model
* `rmodel`  - theoretical model (e.g. `SphericalVariogram`)
* `rcolor`  - color of range curve

## Theoretical variogram options:

* `maxlag` - maximum lag for theoretical model

### Notes

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function varioplot end
function varioplot! end
