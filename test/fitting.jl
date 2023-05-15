@testset "Fitting" begin
  wl = geostatsimage("WalkerLake")
  TI = asarray(wl, :Z)[1:20, 1:20]
  d = georef((z=TI,))
  g = EmpiricalVariogram(d, :z, maxlag=15.0)

  # all fits lead to similar sill
  γ₁ = fit(GaussianVariogram, g)
  γ₂ = fit(SphericalVariogram, g)
  γ₃ = fit(ExponentialVariogram, g)
  γ₄ = fit(MaternVariogram, g)
  @test isapprox(sill(γ₁), 0.054, atol=1e-3)
  @test isapprox(sill(γ₂), 0.054, atol=1e-3)
  @test isapprox(sill(γ₃), 0.058, atol=1e-3)
  @test isapprox(sill(γ₄), 0.057, atol=1e-3)

  # best fit is a Gaussian variogram
  γ = fit(Variogram, g)
  @test γ isa GaussianVariogram
  @test isapprox(sill(γ), 0.054, atol=1e-3)

  # make sure convenient methods work
  γ₁ = fit(GaussianVariogram, g, h -> 1 / h)
  γ₂ = fit(Variogram, g, h -> 1 / h)
  @test sill(γ₁) > 0
  @test sill(γ₂) > 0

  # unitful types
  wl = geostatsimage("WalkerLake")
  TI = asarray(wl, :Z)[1:20, 1:20]
  d = georef((z=TI * u"K",))
  g = EmpiricalVariogram(d, :z, maxlag=15.0)
  γ = fit(Variogram, g)
  @test unit(sill(γ)) == u"K^2"
  @test unit(nugget(γ)) == u"K^2"
end
