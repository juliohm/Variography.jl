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
    @test all(γ(h) .≤ γ(h.+1))
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

  if ismaintainer || istravis
    @testset "Plot recipe" begin
      function plot_variograms(fname)
        plt1 = plot()
        for γ ∈ γs
          plot!(plt1, γ, maxlag=3.)
        end
        plt2 = plot()
        for γ ∈ γn
          plot!(plt2, γ, maxlag=3.)
        end
        plot(plt1, plt2, size=(600,800), layout=(2,1))
        png(fname)
      end
      refimg = joinpath(datadir,"TheoreticalVariograms.png")
      @test test_images(VisualTest(plot_variograms, refimg), popup=!istravis, tol=0.1) |> success
    end
  end
end
