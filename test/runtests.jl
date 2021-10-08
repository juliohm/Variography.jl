using Variography
using Meshes
using GeoStatsBase
using GeoStatsImages
using Distances
using LinearAlgebra
using CSV, DelimitedFiles
using Plots; gr(size=(600,400))
using ReferenceTests, ImageIO
using Test, Random

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
