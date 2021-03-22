@testset "Nested" begin
  # nested variogram with nugget effect
  γ = NuggetEffect(0.2) + GaussianVariogram(nugget=0.1, sill=0.8, range=50.0)
  @test nugget(γ) ≈ 0.3
  @test sill(γ) ≈ 1.0
  @test range(γ) ≈ 50.0
  γ = 2.0*NuggetEffect(0.2)
  @test nugget(γ) ≈ 0.4
  @test sill(γ) ≈ 0.4
  @test range(γ) ≈ 0.0

  # sill is defined for nested models
  γ = GaussianVariogram(sill=1.) + ExponentialVariogram(sill=2.)
  @test sill(γ) == 3.0

  # nugget is defined for nested models
  γ₁ = GaussianVariogram()
  γ₂ = GaussianVariogram() + ExponentialVariogram()
  @test nugget(γ₁) == nugget(γ₂)

  # stationarity of nested models
  γ = GaussianVariogram() + ExponentialVariogram() + SphericalVariogram()
  @test isstationary(γ)
  @test sill(γ) == 3.0
  @test !isstationary(γ + PowerVariogram())

  # result type is defined for nested models
  # see https://github.com/JuliaEarth/GeoStats.jl/issues/121 
  γ = GaussianVariogram() + ExponentialVariogram()
  @test Variography.result_type(γ, rand(Point3), rand(Point3)) == Float64
  γ = GaussianVariogram(sill=1.0f0,range=1.0f0,nugget=0.1f0)
  @test Variography.result_type(γ, rand(Point3f), rand(Point3f)) == Float32

  # nested model with matrix coefficients
  C₁ = [1.0 0.5; 0.5 2.0]
  C₂ = [3.0 0.0; 0.0 3.0]
  γ = C₁ * GaussianVariogram(range=1.0) + C₂ * SphericalVariogram(range=2.0)
  @test range(γ) ≈ 2.0
  @test sill(γ)  ≈ C₁ .+ C₂
  @test γ(10.0) ≈ sill(γ)
  @test γ(Point(10.,0.), Point(0.,0.)) ≈ sill(γ)
  @test isstationary(γ)

  # nested model with matrix coefficients
  C = [1. 0.; 0. 1.]
  γ = C*GaussianVariogram() + C*ExponentialVariogram() + C*CubicVariogram()
  @test range(γ) ≈ 1.0
  @test sill(γ) ≈ [3. 0.; 0. 3.]
  @test γ(10.0) ≈ sill(γ)
  @test γ(Point(10.,0.), Point(0.,0.)) ≈ sill(γ)
  @test isstationary(γ)

  # test constructor explicitly
  γ = NestedVariogram((1., 2.), (ExponentialVariogram(), SphericalVariogram()))
  @test sill(γ) == 3.0
  @test range(γ) == 1.0
  @test nugget(γ) == 0.0

  # test individual structures
  γ = SphericalVariogram() + 2ExponentialVariogram() + NuggetEffect(10.0)
  @test structures(γ) == (10.0, (1.0, 2.0), (SphericalVariogram(), ExponentialVariogram()))
  γ = SphericalVariogram(sill=2.0) + ExponentialVariogram(nugget=0.1)
  @test structures(γ) == (0.1, (2.0, 0.9), (SphericalVariogram(), ExponentialVariogram()))
  @test structures(SphericalVariogram()) == (0.0, (1.0,), (SphericalVariogram(),))
end
