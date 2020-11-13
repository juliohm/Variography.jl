# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(γ::EmpiricalVariogram; showbins=true)
  # get the data
  x, y, n = values(γ)

  binsize = x[2] - x[1]

  # discard empty bins
  x = x[n .> 0]; y = y[n .> 0]; n = n[n .> 0]

  # draw bin counts as a measure of confidence
  if showbins
    @series begin
      # plot a "frequency" instead of raw counts
      f = n*(maximum(y) / maximum(n)) / 10.

      seriestype := :bar
      seriescolor := :black
      fillalpha := .1
      label := "bin counts"

      x, f
    end
  end

  seriestype --> :scatterpath
  linecolor --> :black
  xlims --> (0, maximum(x) + binsize/2)
  xguide --> "h"
  yguide --> "γ(h)"
  label --> "variogram"

  x, y
end
