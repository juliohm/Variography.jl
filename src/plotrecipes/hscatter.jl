# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

@userplot HScatter

@recipe function f(hs::HScatter; lag=0, tol=1e-1, distance=Euclidean())
  # get inputs
  sdata = hs.args[1]
  var₁ = hs.args[2]
  var₂ = length(hs.args) == 3 ? hs.args[3] : var₁

  # lookup valid data
  locs₁ = findall(!ismissing, sdata[var₁])
  locs₂ = findall(!ismissing, sdata[var₂])
  𝒟₁ = view(sdata, locs₁)
  𝒟₂ = view(sdata, locs₂)
  X₁, z₁ = coordinates(𝒟₁), 𝒟₁[var₁]
  X₂, z₂ = coordinates(𝒟₂), 𝒟₂[var₂]

  # compute pairwise distance
  m, n = length(z₁), length(z₂)
  pairs = [(i,j) for j in 1:n for i in j:m]
  ds = [evaluate(distance, view(X₁,:,i), view(X₂,:,j)) for (i,j) in pairs]

  xguide --> var₁
  yguide --> var₂
  legend --> false
  aspect_ratio --> :equal

  # find indices with given lag
  match = findall(abs.(ds .- lag) .< tol)

  if isempty(match)
    @warn "no points were found with lag = $lag, skipping..."
    return nothing
  end

  # scatter plot coordinates
  mpairs = view(pairs, match)
  x = z₁[first.(mpairs)]
  y = z₂[last.(mpairs)]

  # plot identity line
  @series begin
    seriestype := :path
    seriescolor := :black
    primary := false
    linestyle := :dash

    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    vmin = min(xmin, ymin)
    vmax = max(xmax, ymax)

    [vmin, vmax], [vmin, vmax]
  end

  # plot h-scatter
  @series begin
    seriestype := :scatter

    x, y
  end
end
