using Variography
using Meshes
using GeoStatsBase
using GeoStatsImages
using Unitful
using Distances
using LinearAlgebra
using CSV, DelimitedFiles
using Plots; gr(size=(600,400))
using MeshPlots # TODO: replace by MeshViz
using ReferenceTests, ImageIO
using Test, Random
import CoDa: Composition

# workaround for GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__,"data")

# list of tests
testfiles = [
  "empirical.jl",
  "partition.jl",
  "varioplane.jl",
  "theoretical.jl",
  "nesting.jl",
  "sampling.jl",
  "pairwise.jl",
  "fitting.jl"
]

@testset "Variography.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
