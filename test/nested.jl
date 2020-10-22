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
  @test Variography.result_type(γ, rand(3), rand(3)) == Float64
  γ = GaussianVariogram(sill=1.0f0,range=1.0f0,nugget=0.1f0)
  @test Variography.result_type(γ, rand(Float32, 3), rand(Float32, 3)) == Float32

  # nested model with matrix coefficients
  C₁ = [1.0 0.5; 0.5 2.0]
  C₂ = [3.0 0.0; 0.0 3.0]
  γ = C₁ * GaussianVariogram(range=1.0) + C₂ * SphericalVariogram(range=2.0)
  @test range(γ) ≈ 2.0
  @test sill(γ)  ≈ C₁ .+ C₂
  @test γ(10.0) ≈ sill(γ)
  @test γ([10.,0.], [0.,0.]) ≈ sill(γ)
  @test isstationary(γ)

  # nested model with mixed coefficients
  γ = GaussianVariogram() + [1. 0.; 0. 1.] * ExponentialVariogram() + CubicVariogram()
  @test range(γ) ≈ 1.0
  @test sill(γ) ≈ [3. 0.; 0. 3.]
  @test γ(10.0) ≈ sill(γ)
  @test γ([10.,0.], [0.,0.]) ≈ sill(γ)
  @test isstationary(γ)

  # test constructor explicitly
  γ = NestedVariogram((1., 2.), (ExponentialVariogram(), SphericalVariogram()))
  @test sill(γ) == 3.0
  @test range(γ) == 1.0
  @test nugget(γ) == 0.0
end
