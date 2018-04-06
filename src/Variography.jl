# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

__precompile__()

module Variography

importall GeoStatsBase

using Distances
using SpecialFunctions: besselk
using RecipesBase

# won't be neeeded in Julia v0.7
using Parameters

# extend pairwise for theoretical variograms
import Distances: pairwise

# variogram models
include("empirical_variograms.jl")
include("theoretical_variograms.jl")

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
  pairwise

end
