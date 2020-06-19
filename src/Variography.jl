# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Variography

using GeoStatsBase

using Optim
using Distances
using InteractiveUtils: subtypes
using SpecialFunctions: gamma, besselk
using Transducers: Map, reduce
using StaticArrays: MVector
using RecipesBase
using Parameters
using Printf

import Base: merge, +
import GeoStatsBase: fit

include("utils.jl")
include("empirical.jl")
include("partition.jl")
include("varioplane.jl")
include("theoretical.jl")
include("pairwise.jl")
include("fitting.jl")

# plot recipes
include("plotrecipes/empirical.jl")
include("plotrecipes/varioplane.jl")
include("plotrecipes/theoretical.jl")
include("plotrecipes/hscatter.jl")

export
  # empirical variograms
  EmpiricalVariogram,
  DirectionalVariogram,
  EmpiricalVarioplane,

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
  VariogramFitAlgo,
  WeightedLeastSquares

end
