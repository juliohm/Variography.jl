# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Variography

using Meshes
using GeoStatsBase

using Optim
using Tables
using Distances
using Bessels: gamma, besselk
using InteractiveUtils: subtypes
using NearestNeighbors: MinkowskiMetric
using Transducers: Map, foldxt
using LinearAlgebra
using Setfield
using Unitful
using Random
using Printf

import Base: merge, +, *
import GeoStatsBase: fit
import Meshes: isisotropic
import LinearAlgebra: ⋅

include("utils.jl")
include("empirical.jl")
include("partition.jl")
include("varioplane.jl")
include("theoretical.jl")
include("nesting.jl")
include("sampling.jl")
include("pairwise.jl")
include("fitting.jl")

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
  CircularVariogram,
  sill,
  nugget,
  isstationary,
  isisotropic,
  structures,
  pairwise,
  pairwise!,

  # fitting methods
  VariogramFitAlgo,
  WeightedLeastSquares

end
