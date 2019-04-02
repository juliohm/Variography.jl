@testset "Empirical variograms" begin
  # homogeneous field has zero variogram
  γ = EmpiricalVariogram(Matrix(1.0I, 3, 3), ones(3), nlags=2, maxlag=2.)
  x, y, n = values(γ)
  @test x ≈ [1/2, 3/2]
  @test isnan(y[1]) && y[2] == 0.
  @test n == [0, 3]

  # test spatial data interface
  γ = EmpiricalVariogram(psetdata2D, :z, nlags=20, maxlag=1.)
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

  # empirical variogram with only missing data
  X = rand(3,2); z = [NaN, NaN, NaN]
  γ = EmpiricalVariogram(X, z, maxlag=1., nlags=5)
  x, y, n = values(γ)
  @test x == [.1, .3, .5, .7, .9]
  @test all(isnan.(y))
  @test all(iszero.(n))

  # directional variogram and known anisotropy ratio
  imgdata = readdlm(joinpath(datadir,"anisotropic.tsv"))
  geodata = RegularGridData{Float64}(Dict(:z => imgdata))
  γhor = DirectionalVariogram(geodata, (1.,0.), :z, maxlag=50.)
  γver = DirectionalVariogram(geodata, (0.,1.), :z, maxlag=50.)
  γₕ = Variography.fit(GaussianVariogram, γhor)
  γᵥ = Variography.fit(GaussianVariogram, γver)
  @test range(γₕ) / range(γᵥ) ≈ 3. atol=.1

  if visualtests
    @plottest plot(γwalker) joinpath(datadir,"EmpiricalVariograms.png") !istravis

    @plottest begin
      p1 = plot(γhor, showbins=false, label="horizontal")
      plot!(γver, showbins=false, label="vertical")
      p2 = plot(γₕ, maxlag=50., label="horizontal")
      plot!(γᵥ, maxlag=50., label="vertical")
      plot(p1, p2, layout=(2,1))
    end joinpath(datadir,"DirectionalVariograms.png") !istravis
  end
end
