using Variography
using GeoStatsDevTools
using GeoStatsImages
using LinearAlgebra
using Plots; gr(size=(600,400))
# using VisualRegressionTests
using Test

# list of maintainers
maintainers = ["juliohm"]

# environment settings
istravis = "TRAVIS" ∈ keys(ENV)
ismaintainer = "USER" ∈ keys(ENV) && ENV["USER"] ∈ maintainers
datadir = joinpath(@__DIR__,"data")

# load data sets
fname2D = joinpath(datadir,"data2D.tsv")
data2D = readtable(fname2D, delim='\t', coordnames=[:x,:y])

# empirical variograms
TI = training_image("WalkerLake")[1:20,1:20,1]
xwalker = Float64[i for i=1:20 for j=1:20]
ywalker = Float64[j for i=1:20 for j=1:20]
zwalker = Float64[TI[i,j] for i=1:20 for j=1:20]
γwalker = EmpiricalVariogram(hcat(xwalker,ywalker)', zwalker, maxlag=15.)

# list of tests
testfiles = [
  "empirical_variograms.jl",
  "theoretical_variograms.jl",
  "pairwise.jl",
  "fitting.jl"
]

@testset "Variography.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
