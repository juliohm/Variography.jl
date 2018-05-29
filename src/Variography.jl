# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

__precompile__()

module Variography

importall GeoStatsBase

using Optim
using Distances
using SpecialFunctions: besselk
using StaticArrays
using Missings
using RecipesBase

# won't be neeeded in Julia v0.7
using Parameters

# extend methods for theoretical variograms
import Base: +
import Distances: result_type, pairwise

# variogram models
include("empirical_variograms.jl")
include("theoretical_variograms.jl")
include("pairwise.jl")
include("fitting.jl")

# plot recipes
include("plotrecipes/empirical_variograms.jl")
include("plotrecipes/theoretical_variograms.jl")

export
  # empirical variograms
  EmpiricalVariogram,

  # theoretical variograms
  Variogram,
  GaussianVariogram,
  ExponentialVariogram,
  MaternVariogram,
  SphericalVariogram,
  CubicVariogram,
  PentasphericalVariogram,
  PowerVariogram,
  SineHoleVariogram,
  CompositeVariogram,
  isstationary,
  sill,
  pairwise,

  # fitting methods
  WeightedLeastSquares,
  fit

end
