# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Variography

using Meshes
using GeoStatsBase

using Optim
using Distances
using InteractiveUtils: subtypes
using SpecialFunctions: gamma, besselk
using NearestNeighbors: MinkowskiMetric
using Transducers: Map, foldxt
using LinearAlgebra
using RecipesBase
using Parameters
using Setfield
using Random
using Printf

import Base: merge, +, *
import GeoStatsBase: fit

include("utils.jl")
include("empirical.jl")
include("partition.jl")
include("varioplane.jl")
include("theoretical.jl")
include("nesting.jl")
include("sampling.jl")
include("pairwise.jl")
include("fitting.jl")

# plot recipes
include("plotrecipes/empirical.jl")
include("plotrecipes/varioplane.jl")
include("plotrecipes/theoretical.jl")

export
  # empirical variograms
  EmpiricalVariogram,
  EmpiricalVarioplane,
  DirectionalVariogram,
  PlanarVariogram,
  distance,

  # theoretical variograms
  Variogram,
  NuggetEffect,
  GaussianVariogram,
  ExponentialVariogram,
  MaternVariogram,
  SphericalVariogram,
  CubicVariogram,
  PentasphericalVariogram,
  PowerVariogram,
  SineHoleVariogram,
  NestedVariogram,
  isstationary,
  sill, nugget,
  structures,
  distance,
  pairwise,
  pairwise!,

  # fitting methods
  VariogramFitAlgo,
  WeightedLeastSquares

end
