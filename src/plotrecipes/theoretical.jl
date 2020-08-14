# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

_hmax(γ::Variogram) = 3range(γ)
_hmax(γ::PowerVariogram) = 3.

@recipe function f(γ::Variogram, minlag=0, maxlag=_hmax(γ); nlags=100)
  # start at 1e-6 instead of 0 to avoid nugget artifact in plot
  h = range(minlag+1e-6, stop=maxlag, length=nlags)

  seriestype --> :path
  xlims --> (0, maxlag)
  ylims --> (0, Inf)
  xguide --> "Lag h"
  yguide --> "Variogram(h)"
  label --> "variogram"

  h, γ.(h)
end
