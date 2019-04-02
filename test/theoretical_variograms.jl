@testset "Theoretical variograms" begin
  h = range(0, stop=10, length=50)
  x, y = rand(3), rand(3)

  # stationary variogram models
  γs = [GaussianVariogram(), ExponentialVariogram(),
        MaternVariogram(), SphericalVariogram(),
        SphericalVariogram(range=2.), CubicVariogram(),
        PentasphericalVariogram(), SineHoleVariogram()]

  # non-stationary variogram models
  γn = [PowerVariogram(), PowerVariogram(exponent=.4)]

  # non-decreasing variogram models
  γnd = [GaussianVariogram(), ExponentialVariogram(),
         MaternVariogram(), SphericalVariogram(),
         SphericalVariogram(range=2.), CubicVariogram(),
         PentasphericalVariogram(), PowerVariogram()]

  # check stationarity
  @test all(isstationary(γ) for γ ∈ γs)
  @test all(!isstationary(γ) for γ ∈ γn)

  # variograms are symmetric under Euclidean distance
  for γ in (γs ∪ γn ∪ γnd ∪ [sum(γs) + sum(γn) + sum(γnd)])
    @test γ(x, y) ≈ γ(y, x)
  end

  # some variograms are non-decreasing
  for γ in (γnd ∪ [sum(γnd)])
    @test all(γ.(h) .≤ γ.(h.+1))
  end

  # variograms are valid at the origin
  for γ in (γs ∪ γn ∪ γnd)
    @test !isnan(γ(0.)) && !isinf(γ(0.))
  end

  # sill is defined for compositive stationary models
  γ = GaussianVariogram(sill=1.) + ExponentialVariogram(sill=2.)
  @test sill(γ) == 3.

  # composite (additive) models via addition
  γ = GaussianVariogram() + ExponentialVariogram() + SphericalVariogram()
  @test γ isa CompositeVariogram
  @test isstationary(γ)
  @test sill(γ) == 3.
  @test !isstationary(γ + PowerVariogram())

  # ill-conditioned models and nugget regularization
  # see https://github.com/juliohm/GeoStats.jl/issues/29
  X = [93.0 90.0 89.0 94.0 93.0 97.0 95.0 88.0 96.0 98.0
       40.0 33.0 34.0 36.0 30.0 39.0 39.0 28.0 25.0 35.0]

  # ill-conditioned covariance
  γ = GaussianVariogram(range=20.)
  C = sill(γ) .- pairwise(γ, X)
  @test cond(C) > 1000.

  # nugget regularization
  γ = GaussianVariogram(range=20., nugget=0.1)
  C = sill(γ) .- pairwise(γ, X)
  @test γ(0) == 0
  @test γ(1e-6) > 0
  @test cond(C) < 100.

  if visualtests
    @plottest begin
      plt1 = plot()
      for γ ∈ γs
        plot!(plt1, γ, maxlag=3.)
      end
      plt2 = plot()
      for γ ∈ γn
        plot!(plt2, γ, maxlag=3.)
      end
      plot(plt1, plt2, size=(600,800), layout=(2,1))
    end joinpath(datadir,"TheoreticalVariograms.png") !istravis
  end
end
