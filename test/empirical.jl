@testset "Empirical" begin
  @testset "Variogram" begin
    # homogeneous field has zero variogram
    sdata = georef((z=ones(3),), Matrix(1.0I, 3, 3))
    γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2.0)
    x, y, n = values(γ)
    @test x ≈ [1 / 2, √2]
    @test y[2] == 0.0
    @test n == [0, 3]

    # basic test on number of lags
    sdata = georef((z=[1.0, 0.0, 1.0],), [25.0 50.0 75.0; 25.0 75.0 50.0])
    γ = EmpiricalVariogram(sdata, :z, nlags=20, maxlag=1.0)
    x, y, n = values(γ)
    @test length(x) == 20
    @test length(y) == 20
    @test length(n) == 20

    # empirical variogram on integer coordinates
    sdata = georef((z=ones(3),), Matrix(1I, 3, 3))
    γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2, algorithm=:full)
    x, y, n = values(γ)
    @test x ≈ [1 / 2, √2]
    @test y[2] == 0.0
    @test n == [0, 3]

    # empirical variogram with only missing data
    X = rand(2, 3)
    z = Union{Float64,Missing}[missing, missing, missing]
    𝒟 = georef((z=z,), X)
    γ = EmpiricalVariogram(𝒟, :z, maxlag=1.0, nlags=5)
    x, y, n = values(γ)
    @test x == [0.1, 0.3, 0.5, 0.7, 0.9]
    @test all(iszero.(n))

    # accumulation algorithms give the same result
    Random.seed!(2021)
    sdata = georef((z=rand(1000),), rand(3, 1000))
    γ₁ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algorithm=:full)
    γ₂ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algorithm=:ball)
    @test isequal(values(γ₁), values(γ₂))

    # custom distance is recorded
    sdata = georef((z=rand(1000),), rand(2, 1000))
    γ = EmpiricalVariogram(sdata, :z, distance=Haversine(6371.0), algorithm=:full)
    @test distance(γ) == Haversine(6371.0)

    # print methods
    rng = MersenneTwister(123)
    d = georef((z=rand(rng, 100, 100),))
    γ = EmpiricalVariogram(d, :z)
    @test sprint(show, γ) == "EmpiricalVariogram"
    @test sprint(show, MIME"text/plain"(), γ) ==
          "EmpiricalVariogram\n  abscissa: (0.3535533905932738, 13.84261778461484)\n  ordinate: (0.0, 0.08323850196902784)\n  N° pairs: 2790126"

    # test variography with compositional data
    data = georef((z=rand(Composition{3}, 100),), rand(2, 100))
    γ = EmpiricalVariogram(data, :z, maxlag=1.0, algorithm=:full)
    x, y, n = values(γ)
    @test all(≥(0), x)
    @test all(≥(0), y)
    @test all(>(0), n)

    # test variography with unitful data
    data = georef((z=[1 * u"K" for i in 1:100],), rand(2, 100))
    γ = EmpiricalVariogram(data, :z, nlags=20)
    x, y, n = values(γ)
    @test all(≥(0), x)
    @test y == fill(0.0 * u"K^2", 20)

    # Matheron's vs Cressie's estimator
    data = geostatsimage("Gaussian30x10")
    γ₁ = EmpiricalVariogram(data, :Z, maxlag=50.0, estimator=:matheron)
    γ₂ = EmpiricalVariogram(data, :Z, maxlag=50.0, estimator=:cressie)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test all(isapprox.(y₁, y₂, atol=0.1))
    @test n₁ == n₂

    # specify variables as strings
    data = geostatsimage("Gaussian30x10")
    γ = EmpiricalVariogram(data, "Z", maxlag=50.0)
    x, y, n = values(γ)
    @test all(≥(0), x)
    @test all(>(0.8), y[11:end])
    @test all(≥(0), n)
  end

  @testset "Varioplane" begin
    img = readdlm(joinpath(datadir, "anisotropic.tsv"))
    data = georef((z=img,))
    γ = EmpiricalVarioplane(data, :z, maxlag=50.0)
    @test sprint(show, γ) == "EmpiricalVarioplane"
    @test sprint(show, MIME"text/plain"(), γ) ==
          "EmpiricalVarioplane\n  N° pairs\n  └─0.00° → 372500\n  └─3.67° → 304782\n  └─7.35° → 298306\n  └─11.02° → 297432\n  └─14.69° → 297243\n  ⋮\n  └─165.31° → 293643\n  └─168.98° → 295850\n  └─172.65° → 296931\n  └─176.33° → 306528\n  └─180.00° → 372500"
  end

  @testset "Directional" begin
    # merge operation does not produce NaN
    dir = (0.286788, -0.496732, -0.819152)
    𝒟 = georef(CSV.File(joinpath(datadir, "nanlags.csv")), (:X, :Y, :Z))
    γ = DirectionalVariogram(dir, 𝒟, :Cu, dtol=45, maxlag=150, nlags=20)
    x, y, n = values(γ)
    @test !any(isnan.(x))
    @test !any(isnan.(y))
    @test !any(isnan.(n))

    # directional variogram and known anisotropy ratio
    img = readdlm(joinpath(datadir, "anisotropic.tsv"))
    sdata = georef((z=img,))
    γhor = DirectionalVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    γver = DirectionalVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    γₕ = fit(GaussianVariogram, γhor)
    γᵥ = fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3.0 atol = 0.1
  end

  @testset "Planar" begin
    # directional equals planar rotated by 90 degrees in 2D
    img = readdlm(joinpath(datadir, "anisotropic.tsv"))
    sdata = georef((z=img,))
    γ₁ = PlanarVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    γ₂ = DirectionalVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂
    γ₁ = PlanarVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    γ₂ = DirectionalVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂

    # planar variogram and known anisotropy ratio
    γhor = PlanarVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    γver = PlanarVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    γₕ = fit(GaussianVariogram, γhor)
    γᵥ = fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3.0 atol = 0.1
  end
end
