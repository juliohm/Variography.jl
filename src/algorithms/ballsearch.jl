# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSearchAccum(maxlag, nlags, distance, estimator)

Accumulate pairs of points in geospatial data with
nearest neighbors inside metric ball.
"""
struct BallSearchAccum{T,D,E} <: VariogramAccumAlgo
  maxlag::T
  nlags::Int
  distance::D
  estimator::E
end

function accumulate(data, var₁, var₂, algo::BallSearchAccum)
  # retrieve algorithm parameters
  maxlag = algo.maxlag
  nlags = algo.nlags
  distance = algo.distance
  estimator = algo.estimator

  # retrieve table and point set
  𝒯 = values(data)
  𝒫 = domain(data)

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

  # fast ball search
  ball = MetricBall(maxlag, distance)
  searcher = BallSearch(𝒫, ball)

  # loop over points inside ball
  @inbounds for j in 1:nelements(𝒫)
    pⱼ = 𝒫[j]
    z₁ⱼ = Z₁[j]
    z₂ⱼ = Z₂[j]
    for i in search(pⱼ, searcher)
      i ≤ j && continue # avoid double counting
      pᵢ = 𝒫[i]
      z₁ᵢ = Z₁[i]
      z₂ᵢ = Z₂[i]

      # evaluate geospatial lag
      h = evaluate(distance, coordinates(pᵢ), coordinates(pⱼ))

      # evaluate (cross-)variance
      v = (z₁ᵢ - z₁ⱼ) ⋅ (z₂ᵢ - z₂ⱼ)

      # bin (or lag) where to accumulate result
      lag = ceil(Int, h / δh)
      lag == 0 && @warn "duplicate coordinates found"

      if 0 < lag ≤ nlags && !ismissing(v)
        xsums[lag] += h
        ysums[lag] += v
        counts[lag] += 1
      end
    end
  end

  xsums, ysums, counts
end
