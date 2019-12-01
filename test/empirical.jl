@testset "Empirical" begin
  # homogeneous field has zero variogram
  sdata = PointSetData(Dict(:z=>ones(3)), Matrix(1.0I, 3, 3))
  γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2.)
  x, y, n = values(γ)
  @test x ≈ [1/2, 3/2]
  @test isnan(y[1]) && y[2] == 0.
  @test n == [0, 3]

  # basic test on number of lags
  sdata = PointSetData(Dict(:z=>[1.,0.,1.]), [25. 50. 75.; 25. 75. 50.])
  γ = EmpiricalVariogram(sdata, :z, nlags=20, maxlag=1.)
  x, y, n = values(γ)
  @test length(x) == 20
  @test length(y) == 20
  @test length(n) == 20

  # empirical variogram on integer coordinates
  sdata = PointSetData(Dict(:z=>ones(3)), Matrix(1I, 3, 3))
  γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2)
  x, y, n = values(γ)
  @test x ≈ [1/2, 3/2]
  @test isnan(y[1]) && y[2] == 0.
  @test n == [0, 3]

  # empirical variogram with only missing data
  X = rand(2,3)
  for z in [[NaN,NaN,NaN], [missing,missing,missing]]
    sdata = PointSetData(Dict(:z=>z), X)
    γ = EmpiricalVariogram(sdata, :z, maxlag=1., nlags=5)
    x, y, n = values(γ)
    @test x == [.1, .3, .5, .7, .9]
    @test all(isnan.(y))
    @test all(iszero.(n))
  end

  # accumulation algorithms give the same result
  sdata = PointSetData(Dict(:z=>rand(1000)), rand(3,1000))
  γ₁ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algo=:full)
  γ₂ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algo=:ball)
  @test isequal(values(γ₁), values(γ₂))

  # directional variogram and known anisotropy ratio
  img = readdlm(joinpath(datadir,"anisotropic.tsv"))
  sdata = RegularGridData{Float64}(Dict(:z => img))
  γhor = DirectionalVariogram(sdata, (1.,0.), :z, maxlag=50.)
  γver = DirectionalVariogram(sdata, (0.,1.), :z, maxlag=50.)
  γₕ = fit(GaussianVariogram, γhor)
  γᵥ = fit(GaussianVariogram, γver)
  @test range(γₕ) / range(γᵥ) ≈ 3. atol=.1

  if visualtests
    TI = training_image("WalkerLake")[1:20,1:20,1]
    d = RegularGridData{Float64}(Dict(:z=>TI))
    γ = EmpiricalVariogram(d, :z, maxlag=15.)
    @plottest plot(γ) joinpath(datadir,"EmpiricalVariograms.png") !istravis

    @plottest begin
      p1 = plot(γhor, showbins=false, label="horizontal")
      plot!(γver, showbins=false, label="vertical")
      p2 = plot(γₕ, maxlag=50., label="horizontal")
      plot!(γᵥ, maxlag=50., label="vertical")
      plot(p1, p2, layout=(2,1))
    end joinpath(datadir,"DirectionalVariograms.png") !istravis
  end
end
