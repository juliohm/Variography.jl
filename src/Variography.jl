# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Variography

using GeoStatsBase

using Printf
using Optim
using Distances
using Statistics
using SpecialFunctions: besselk, gamma
using StaticArrays
using RecipesBase
using Parameters

import Base: +
import GeoStatsBase: fit

# variogram models
include("empirical.jl")
include("theoretical.jl")
include("pairwise.jl")
include("fitting.jl")

# plot recipes
include("plotrecipes/empirical.jl")
include("plotrecipes/theoretical.jl")
include("plotrecipes/hscatter.jl")
include("plotrecipes/varplane.jl")

export
  # empirical variograms
  EmpiricalVariogram,
  DirectionalVariogram,

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
  sill, nugget,
  pairwise,
  pairwise!,

  # fitting methods
  WeightedLeastSquares

end
