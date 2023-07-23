# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FullSearchAccum(maxlag, nlags, distance, estimator)

Accumulate pairs of points in geospatial data with
exhaustive (or full) search.
"""
struct FullSearchAccum{T,D,E} <: VariogramAccumAlgo
  maxlag::T
  nlags::Int
  distance::D
  estimator::E
end

neighfun(pset, ::FullSearchAccum) = j -> (j + 1):nelements(pset)

skipfun(::FullSearchAccum) = (i, j) -> false

exitfun(algo::FullSearchAccum) = h -> h > algo.maxlag
