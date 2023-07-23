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

include("algorithms/fullsearch.jl")
include("algorithms/ballsearch.jl")

# ---------------------
# EMPIRICAL VARIOGRAMS
# ---------------------

include("empirical/variogram.jl")
include("empirical/varioplane.jl")
include("empirical/partition.jl")
