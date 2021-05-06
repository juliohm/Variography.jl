@testset "Partition" begin
  @testset "Directional" begin
    # merge operation does not produce NaN
    dir = (0.286788, -0.496732, -0.819152)
    𝒟 = georef(CSV.File(joinpath(datadir,"nanlags.csv")), (:X,:Y,:Z))
    γ = DirectionalVariogram(dir, 𝒟, :Cu, dtol=45, maxlag=150, nlags=20)
    x, y, n = values(γ)
    @test !any(isnan.(x))
    @test !any(isnan.(y))
    @test !any(isnan.(n))

    # directional variogram and known anisotropy ratio
    img = readdlm(joinpath(datadir,"anisotropic.tsv"))
    sdata = georef((z=img,))
    γhor = DirectionalVariogram((1.,0.), sdata, :z, maxlag=50.)
    γver = DirectionalVariogram((0.,1.), sdata, :z, maxlag=50.)
    γₕ = fit(GaussianVariogram, γhor)
    γᵥ = fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3. atol=.1

    if visualtests
      p1 = plot(γhor, showbins=false, label="horizontal")
      plot!(γver, showbins=false, label="vertical")
      p2 = plot(γₕ, 0., 50., label="horizontal")
      plot!(γᵥ, 0., 50., label="vertical")
      plt = plot(p1, p2, layout=(2,1))
      @test_reference "data/directional.png" plt
    end
  end

  @testset "Planar" begin
    # directional equals planar rotated by 90 degrees in 2D
    img = readdlm(joinpath(datadir,"anisotropic.tsv"))
    sdata = georef((z=img,))
    γ₁ = PlanarVariogram((0.,1.), sdata, :z, maxlag=50.)
    γ₂ = DirectionalVariogram((1.,0.), sdata, :z, maxlag=50.)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂
    γ₁ = PlanarVariogram((1.,0.), sdata, :z, maxlag=50.)
    γ₂ = DirectionalVariogram((0.,1.), sdata, :z, maxlag=50.)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂

    # planar variogram and known anisotropy ratio
    γhor = PlanarVariogram((0.,1.), sdata, :z, maxlag=50.)
    γver = PlanarVariogram((1.,0.), sdata, :z, maxlag=50.)
    γₕ = fit(GaussianVariogram, γhor)
    γᵥ = fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3. atol=.1

    if visualtests
      p1 = plot(γhor, showbins=false, label="horizontal")
      plot!(γver, showbins=false, label="vertical")
      p2 = plot(γₕ, 0., 50., label="horizontal")
      plot!(γᵥ, 0., 50., label="vertical")
      plt = plot(p1, p2, layout=(2,1))
      @test_reference "data/planar.png" plt
    end
  end
end
