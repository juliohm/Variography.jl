@testset "Empirical variograms" begin
  # homogeneous field has zero variogram
  γ = EmpiricalVariogram(Matrix(1.0I, 3, 3), ones(3), nlags=2, maxlag=2.)
  x, y, n = values(γ)
  @test x ≈ [1/2, 3/2]
  @test isnan(y[1]) && y[2] == 0.
  @test n == [0, 3]

  # test geodataframe interface
  γ = EmpiricalVariogram(data2D, :value, nlags=20, maxlag=1.)
  x, y, n = values(γ)
  @test length(x) == 20
  @test length(y) == 20
  @test length(n) == 20

  # empirical variogram on integer coordinates
  γ = EmpiricalVariogram(Matrix(1I, 3, 3), ones(3), nlags=2, maxlag=2)
  x, y, n = values(γ)
  @test x ≈ [1/2, 3/2]
  @test isnan(y[1]) && y[2] == 0.
  @test n == [0, 3]

  if ismaintainer || istravis
    @testset "Plot recipe" begin
      function plot_variograms(fname)
        plot(γwalker)
        png(fname)
      end
      refimg = joinpath(datadir,"EmpiricalVariograms.png")
      @test test_images(VisualTest(plot_variograms, refimg), popup=!istravis) |> success
    end
  end
end
