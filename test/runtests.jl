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

# load data sets
fname2D = joinpath(datadir,"data2D.tsv")
data2D = readtable(fname2D, delim='\t', coordnames=[:x,:y])

# empirical variograms
TI = training_image("WalkerLake")[1:20,1:20,1]
x = Float64[i for i=1:20 for j=1:20]
y = Float64[j for i=1:20 for j=1:20]
z = Float64[TI[i,j] for i=1:20 for j=1:20]
γwalker = EmpiricalVariogram(hcat(x,y)', z, maxlag=15.)

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
