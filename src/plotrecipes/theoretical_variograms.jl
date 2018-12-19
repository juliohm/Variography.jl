# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

@recipe function f(γ::Variogram; maxlag=3.)
  # start at 1e-6 instead of 0 to avoid nugget artifact in plot
  h = range(1e-6, stop=maxlag, length=100)

  seriestype --> :path
  xlim --> (0, maxlag)
  xlabel --> "Lag h"
  ylabel --> "Variogram(h)"
  label --> "variogram"

  h, γ(h)
end
