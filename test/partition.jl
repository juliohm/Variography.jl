@testset "Partition" begin
  @testset "Directional" begin
    # directional variogram and known anisotropy ratio
    img = readdlm(joinpath(datadir,"anisotropic.tsv"))
    sdata = RegularGridData(OrderedDict(:z => img))
    γhor = DirectionalVariogram(sdata, (1.,0.), :z, maxlag=50.)
    γver = DirectionalVariogram(sdata, (0.,1.), :z, maxlag=50.)
    γₕ = fit(GaussianVariogram, γhor)
    γᵥ = fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3. atol=.1

    if visualtests
      @plottest begin
        p1 = plot(γhor, showbins=false, label="horizontal")
        plot!(γver, showbins=false, label="vertical")
        p2 = plot(γₕ, maxlag=50., label="horizontal")
        plot!(γᵥ, maxlag=50., label="vertical")
        plot(p1, p2, layout=(2,1))
      end joinpath(datadir,"directional.png") !istravis
    end
  end

  @testset "Planar" begin
    # directional equals planar rotated by 90 degrees in 2D
    img = readdlm(joinpath(datadir,"anisotropic.tsv"))
    sdata = RegularGridData(OrderedDict(:z => img))
    γ₁ = PlanarVariogram(sdata, (0.,1.), :z, maxlag=50.)
    γ₂ = DirectionalVariogram(sdata, (1.,0.), :z, maxlag=50.)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂
    γ₁ = PlanarVariogram(sdata, (1.,0.), :z, maxlag=50.)
    γ₂ = DirectionalVariogram(sdata, (0.,1.), :z, maxlag=50.)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂

    # planar variogram and known anisotropy ratio
    γhor = PlanarVariogram(sdata, (0.,1.), :z, maxlag=50.)
    γver = PlanarVariogram(sdata, (1.,0.), :z, maxlag=50.)
    γₕ = fit(GaussianVariogram, γhor)
    γᵥ = fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3. atol=.1

    if visualtests
      @plottest begin
        p1 = plot(γhor, showbins=false, label="horizontal")
        plot!(γver, showbins=false, label="vertical")
        p2 = plot(γₕ, maxlag=50., label="horizontal")
        plot!(γᵥ, maxlag=50., label="vertical")
        plot(p1, p2, layout=(2,1))
      end joinpath(datadir,"planar.png") !istravis
    end
  end
end
