# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

module Variography

using GeoStatsBase
using GeoStatsDevTools

using Printf
using Optim
using Distances
using Statistics
using SpecialFunctions: besselk, gamma
using StaticArrays
using RecipesBase
using Parameters

# extend methods for theoretical variograms
import Base: +

# variogram models
include("empirical_variograms.jl")
include("theoretical_variograms.jl")
include("pairwise.jl")
include("fitting.jl")

# plot recipes
include("plotrecipes/empirical_variograms.jl")
include("plotrecipes/theoretical_variograms.jl")
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
  sill,
  pairwise,
  pairwise!,

  # fitting methods
  WeightedLeastSquares,
  fit

end
