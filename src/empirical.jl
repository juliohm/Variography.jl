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
    AccumAlgorithm

An accumulation algorithm for estimating [`EmpiricalVariogram`](@ref).
"""
abstract type VariogramAccumAlgo end

"""
    accumulate(data, var₁, var₂, algo)

Accumulate pairs of points in `data` for variables
`var₁` and `var₂` with accumulation algorithm `algo`.
"""
function accumulate end

include("algorithms/fullsearch.jl")
include("algorithms/ballsearch.jl")

# ---------------------
# EMPIRICAL VARIOGRAMS
# ---------------------

include("empirical/variogram.jl")
include("empirical/varioplane.jl")
include("empirical/partition.jl")
