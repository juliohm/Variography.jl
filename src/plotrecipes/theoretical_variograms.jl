# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

@recipe function f(γ::Variogram; nlags=100, maxlag=3.)
  # start at 1e-6 instead of 0 to avoid nugget artifact in plot
  h = range(1e-6, stop=maxlag, length=nlags)

  seriestype --> :path
  xlim --> (0, maxlag)
  xguide --> "Lag h"
  yguide --> "Variogram(h)"
  label --> "variogram"

  h, γ.(h)
end
