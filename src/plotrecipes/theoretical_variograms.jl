# ------------------------------------------------------------------
# Copyright (c) 2017, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

@recipe function f(γ::Variogram; maxlag=3.)
  # discretize
  h = range(0, stop=maxlag, length=100)

  seriestype --> :path
  xlabel --> "Lag h"
  ylabel --> "Variogram(h)"
  label --> "variogram"

  h, γ(h)
end
