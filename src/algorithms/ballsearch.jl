# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSearchAccum(maxlag, nlags, distance)

Accumulate pairs of points in geospatial data with
nearest neighbors inside metric ball.
"""
struct BallSearchAccum{T,D} <: VariogramAccumAlgo
  nlags::Int
  maxlag::T
  distance::D
end

function _neighbors(algo::BallSearchAccum, pset, j)
  ball = MetricBall(algo.maxlag, algo.distance)
  searcher = BallSearch(pset, ball)
  @inbounds(search(pset[j], searcher))
end

_skip(::BallSearchAccum, i, j) = i ≤ j

_exit(::BallSearchAccum, h) = false
