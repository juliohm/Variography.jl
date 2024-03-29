@testset "Fitting" begin
  img = readdlm(joinpath(datadir, "WalkerLake.txt"))
  d = georef((; Z=img))
  g = EmpiricalVariogram(d, :Z, maxlag=15.0)

  # all fits lead to similar sill
  γ₁ = Variography.fit(GaussianVariogram, g)
  γ₂ = Variography.fit(SphericalVariogram, g)
  γ₃ = Variography.fit(ExponentialVariogram, g)
  γ₄ = Variography.fit(MaternVariogram, g)
  @test isapprox(sill(γ₁), 0.054, atol=1e-3)
  @test isapprox(sill(γ₂), 0.054, atol=1e-3)
  @test isapprox(sill(γ₃), 0.054, atol=1e-3)
  @test isapprox(sill(γ₄), 0.054, atol=1e-3)

  # fix parameters
  γ = Variography.fit(GaussianVariogram, g, range=12.0)
  @test isapprox(range(γ), 12.0, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, sill=0.07)
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, nugget=0.05)
  @test isapprox(nugget(γ), 0.05, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, range=12.0, sill=0.07)
  @test isapprox(range(γ), 12.0, atol=1e-3)
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, range=12.0, nugget=0.05)
  @test isapprox(range(γ), 12.0, atol=1e-3)
  @test isapprox(nugget(γ), 0.05, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, sill=0.07, nugget=0.05)
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  @test isapprox(nugget(γ), 0.05, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, range=12.0, sill=0.07, nugget=0.05)
  @test isapprox(range(γ), 12.0, atol=1e-3)
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  @test isapprox(nugget(γ), 0.05, atol=1e-3)

  # fix maximum parameters
  γ = Variography.fit(GaussianVariogram, g, maxrange=5.0)
  @test isapprox(range(γ), 5.0, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, maxsill=0.04)
  @test isapprox(sill(γ), 0.04, atol=1e-3)
  γ = Variography.fit(GaussianVariogram, g, maxnugget=0.004)
  @test isapprox(nugget(γ), 0.004, atol=1e-3)

  # best fit is a Gaussian variogram
  γ = Variography.fit(Variogram, g)
  @test γ isa GaussianVariogram
  @test isapprox(sill(γ), 0.054, atol=1e-3)
  γ = Variography.fit([SphericalVariogram, GaussianVariogram], g)
  @test γ isa GaussianVariogram
  @test isapprox(sill(γ), 0.054, atol=1e-3)

  # make sure convenient methods work
  γ₁ = Variography.fit(GaussianVariogram, g, h -> 1 / h)
  γ₂ = Variography.fit(Variogram, g, h -> 1 / h)
  @test sill(γ₁) > 0
  @test sill(γ₂) > 0

  # unitful types
  img = readdlm(joinpath(datadir, "WalkerLake.txt"))
  d = georef((; Z=img * u"K"))
  g = EmpiricalVariogram(d, :Z, maxlag=15.0)
  γ = Variography.fit(Variogram, g)
  @test unit(sill(γ)) == u"K^2"
  @test unit(nugget(γ)) == u"K^2"
end
