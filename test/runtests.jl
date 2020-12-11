using Variography
using GeoStatsBase
using GeoStatsImages
using Distances
using LinearAlgebra
using DelimitedFiles
using Plots; gr(size=(600,400))
using VisualRegressionTests
using Test, Pkg, Random

# workaround for GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
if !isCI
  Pkg.add("Gtk")
  using Gtk
end
datadir = joinpath(@__DIR__,"data")

# list of tests
testfiles = [
  "empirical.jl",
  "partition.jl",
  "varioplane.jl",
  "theoretical.jl",
  "nested.jl",
  "pairwise.jl",
  "fitting.jl"
]

@testset "Variography.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
