using Variography
using GeoStatsBase
using GeoStatsImages
using Distances
using LinearAlgebra
using DelimitedFiles
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

# helper functions for visual regression tests
function asimage(plt)
  io = IOBuffer()
  show(io, "image/png", plt)
  seekstart(io)
  ImageIO.load(io)
end
macro test_ref_plot(fname, plt)
  esc(quote
    @test_reference $fname asimage($plt)
  end)
end

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
