@testset "Fitting" begin
  # variogram types to fit
  Vs = (GaussianVariogram, SphericalVariogram,
        ExponentialVariogram, MaternVariogram)

  # all fits lead to same sill
  γs = [fit(V, γwalker) for V in Vs]
  for γ in γs
    @test isapprox(sill(γ), 0.054, atol=1e-3)
  end

  # best fit is a Gaussian variogram
  γbest = fit(Variogram, γwalker)
  @test γbest isa GaussianVariogram
  @test isapprox(sill(γbest), 0.054, atol=1e-3)

  if ismaintainer || istravis
    @testset "Plot recipe" begin
      function plot_fit(fname)
        plts = []
        for γ in γs
          plt = plot(γwalker, legend=false)
          plot!(γ, maxlag=15.)
          push!(plts, plt)
        end
        plot(plts...)
        png(fname)
      end
      refimg = joinpath(datadir,"Fitting.png")
      @test test_images(VisualTest(plot_fit, refimg), popup=!istravis) |> success
    end
  end
end
