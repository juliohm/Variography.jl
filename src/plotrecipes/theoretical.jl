# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

_hmax(γ::Variogram) = 3range(γ)
_hmax(γ::PowerVariogram) = 3.0
_hmax(γ::NuggetEffect) = 3.0

@recipe function f(γ::Variogram, minlag=0, maxlag=_hmax(γ); nlags=100)
  # start at 1e-6 instead of 0 to avoid nugget artifact in plot
  h = range(minlag+1e-6, stop=maxlag, length=nlags)

  seriestype --> :path
  label --> string(nameof(typeof(γ)))
  xlims --> (0, maxlag)
  ylims --> (0, Inf)
  xguide --> "h"
  yguide --> "γ(h)"

  h, γ.(h)
end
