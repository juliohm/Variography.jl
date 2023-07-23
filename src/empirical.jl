# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    VariogramEstimator

A (robust) estimator of [`EmpiricalVariogram`](@ref).
"""
abstract type VariogramEstimator end

result_type(estim::VariogramEstimator, z₁, z₂) = typeof(formula(estim, z₁[1], z₁[2], z₂[1], z₂[2]))

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
function accumulate(data, var₁, var₂, estim::VariogramEstimator, algo::VariogramAccumAlgo)
  # retrieve algorithm parameters
  nlags = algo.nlags
  maxlag = algo.maxlag
  distance = algo.distance

  # compute lag size
  δh = maxlag / nlags

  # table and point set
  𝒯 = values(data)
  𝒫 = domain(data)

  # vectors for variables
  cols = Tables.columns(𝒯)
  z₁ = Tables.getcolumn(cols, var₁)
  z₂ = Tables.getcolumn(cols, var₂)

  # neighbors function
  neighbors = neighfun(algo, 𝒫)

  # skip condition
  skip = skipfun(algo)

  # early exit condition
  exit = exitfun(algo)

  # accumulation type
  V = result_type(estim, z₁, z₂)

  # lag sums and counts
  xsums = zeros(nlags)
  ysums = zeros(V, nlags)
  counts = zeros(Int, nlags)

  # loop over points inside ball
  @inbounds for j in 1:nelements(𝒫)
    pⱼ = 𝒫[j]
    z₁ⱼ = z₁[j]
    z₂ⱼ = z₂[j]
    for i in neighbors(j)
      # skip to avoid double counting
      skip(i, j) && continue

      pᵢ = 𝒫[i]
      z₁ᵢ = z₁[i]
      z₂ᵢ = z₂[i]

      # evaluate geospatial lag
      h = evaluate(distance, coordinates(pᵢ), coordinates(pⱼ))

      # early exit if out of range
      exit(h) && continue

      # evaluate (cross-)variance
      v = formula(estim, z₁ᵢ, z₁ⱼ, z₂ᵢ, z₂ⱼ)

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

  # bin (or lag) size
  lags = range(δh / 2, stop=maxlag - δh / 2, length=nlags)

  # ordinate function
  ordfun(ysum, count) = normsum(estim, ysum, count)

  # variogram abscissa
  abscissa = @. xsums / counts
  abscissa[counts .== 0] .= lags[counts .== 0]

  # variogram ordinate
  ordinate = @. ordfun(ysums, counts)
  ordinate[counts .== 0] .= zero(eltype(ordinate))

  abscissa, ordinate, counts
end

include("algorithms/fullsearch.jl")
include("algorithms/ballsearch.jl")

# ---------------------
# EMPIRICAL VARIOGRAMS
# ---------------------

include("empirical/variogram.jl")
include("empirical/varioplane.jl")
include("empirical/partition.jl")
