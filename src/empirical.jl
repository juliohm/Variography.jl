# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    VariogramEstimator

A (robust) estimator of [`EmpiricalVariogram`](@ref).
"""
abstract type VariogramEstimator end

include("estimators/matheron.jl")

"""
    VariogramAccumAlgo

Algorithm for accumulating pairs of points in
[`EmpiricalVariogram`](@ref) estimation.
"""
abstract type VariogramAccumAlgo end

"""
    accumulate(data, var₁, var₂, algo)

Accumulate pairs of points in `data` for variables
`var₁` and `var₂` with accumulation algorithm `algo`.
"""
function accumulate(data, var₁, var₂, algo::VariogramAccumAlgo)
  # retrieve algorithm parameters
  maxlag = algo.maxlag
  nlags = algo.nlags
  distance = algo.distance
  estimator = algo.estimator

  # retrieve table and point set
  𝒯 = values(data)
  𝒫 = domain(data)

  # retrieve neighbors function
  neighbors = neighfun(𝒫, algo)

  # retrieve skip condition
  skip = skipfun(algo)

  # retrieve early exit condition
  exit = exitfun(algo)

  # collect vectors for variables
  cols = Tables.columns(𝒯)
  Z₁ = Tables.getcolumn(cols, var₁)
  Z₂ = Tables.getcolumn(cols, var₂)

  # lag size
  δh = maxlag / nlags

  # accumulation type
  V = typeof((Z₁[1] - Z₂[1]) ⋅ (Z₁[1] - Z₂[1]))

  # lag sums and counts
  xsums = zeros(nlags)
  ysums = zeros(V, nlags)
  counts = zeros(Int, nlags)

  # loop over points inside ball
  @inbounds for j in 1:nelements(𝒫)
    pⱼ = 𝒫[j]
    z₁ⱼ = Z₁[j]
    z₂ⱼ = Z₂[j]
    for i in neighbors(j)
      # skip to avoid double counting
      skip(i, j) && continue

      pᵢ = 𝒫[i]
      z₁ᵢ = Z₁[i]
      z₂ᵢ = Z₂[i]

      # evaluate geospatial lag
      h = evaluate(distance, coordinates(pᵢ), coordinates(pⱼ))

      # early exit if out of range
      exit(h) && continue

      # evaluate (cross-)variance
      v = (z₁ᵢ - z₁ⱼ) ⋅ (z₂ᵢ - z₂ⱼ)

      # bin (or lag) where to accumulate result
      lag = ceil(Int, h / δh)
      lag == 0 && @warn "duplicate coordinates found, consider using `UniqueCoords`"

      if 0 < lag ≤ nlags && !ismissing(v)
        xsums[lag] += h
        ysums[lag] += v
        counts[lag] += 1
      end
    end
  end

  xsums, ysums, counts
end

include("algorithms/fullsearch.jl")
include("algorithms/ballsearch.jl")

# ---------------------
# EMPIRICAL VARIOGRAMS
# ---------------------

include("empirical/variogram.jl")
include("empirical/varioplane.jl")
include("empirical/partition.jl")
