@testset "Partition" begin
  @testset "Directional" begin
    # directional variogram and known anisotropy ratio
    img = readdlm(joinpath(datadir,"anisotropic.tsv"))
    sdata = RegularGridData{Float64}(OrderedDict(:z => img))
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
end
