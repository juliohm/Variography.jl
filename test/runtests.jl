using Variography
using Meshes
using GeoTables
using Unitful
using Distances
using LinearAlgebra
using CSV, DelimitedFiles
using Test, Random
import CoDa: Composition

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = ["empirical.jl", "theoretical.jl", "nesting.jl", "fitting.jl", "sampling.jl"]

@testset "Variography.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
