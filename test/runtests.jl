using Variography
using Meshes
using GeoStatsBase
using GeoStatsImages
using Unitful
using Distances
using LinearAlgebra
using CSV, DelimitedFiles
using Test, Random
import CoDa: Composition

# environment settings
datadir = joinpath(@__DIR__, "data")

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
