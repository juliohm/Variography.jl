using Variography
using GeoStatsDevTools
using GeoStatsImages
using Plots; gr(size=(600,400))
using Base.Test
using VisualRegressionTests

# setup GR backend for Travis CI
ENV["GKSwstype"] = "100"
ENV["PLOTS_TEST"] = "true"

# list of maintainers
maintainers = ["juliohm"]

# environment settings
istravis = "TRAVIS" ∈ keys(ENV)
ismaintainer = "USER" ∈ keys(ENV) && ENV["USER"] ∈ maintainers
datadir = joinpath(@__DIR__,"data")

fname2D = joinpath(datadir,"data2D.tsv")
data2D = readtable(fname2D, delim='\t', coordnames=[:x,:y])

# list of tests
testfiles = [
  "empirical_variograms.jl",
  "theoretical_variograms.jl"
]

@testset "Variography.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
