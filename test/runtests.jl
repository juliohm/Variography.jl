using Variography
using GeoStatsBase
using GeoStatsImages
using LinearAlgebra
using DelimitedFiles
using Plots; gr(size=(600,400))
using VisualRegressionTests
using Test, Pkg

# workaround for GR warnings
ENV["GKSwstype"] = "100"

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
datadir = joinpath(@__DIR__,"data")
visualtests = !istravis || (istravis && islinux)
if !istravis
  Pkg.add("Gtk")
  using Gtk
end

# list of tests
testfiles = [
  "empirical_variograms.jl",
  "theoretical_variograms.jl",
  "pairwise.jl",
  "fitting.jl",
  "plotrecipes.jl"
]

@testset "Variography.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
