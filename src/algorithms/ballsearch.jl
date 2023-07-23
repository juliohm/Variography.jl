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

function neighfun(pset, algo::BallSearchAccum)
  ball = MetricBall(algo.maxlag, algo.distance)
  searcher = BallSearch(pset, ball)
  j -> @inbounds(search(pset[j], searcher))
end

skipfun(::BallSearchAccum) = (i, j) -> i ≤ j

exitfun(::BallSearchAccum) = h -> false
