# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FullSearchAccum(maxlag, nlags, distance)

Accumulate pairs of points in geospatial data with
exhaustive (or full) search.
"""
struct FullSearchAccum{T,D} <: VariogramAccumAlgo
  nlags::Int
  maxlag::T
  distance::D
end

neighfun(::FullSearchAccum, pset) = j -> (j + 1):nelements(pset)

skipfun(::FullSearchAccum) = (i, j) -> false

exitfun(algo::FullSearchAccum) = h -> h > algo.maxlag
