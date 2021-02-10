@testset "Varioplane" begin
  Random.seed!(123)
  img = readdlm(joinpath(datadir,"anisotropic.tsv"))
  sdata = georef((z=img,))
  γ = EmpiricalVarioplane(sdata, :z, maxlag=50.)

  @test sprint(show, γ) == "EmpiricalVarioplane"
  @test sprint(show, MIME"text/plain"(), γ) == "EmpiricalVarioplane\n  N° pairs\n  └─0.00° → 372500\n  └─3.67° → 304782\n  └─7.35° → 298306\n  └─11.02° → 297432\n  └─14.69° → 297243\n  ⋮\n  └─165.31° → 293643\n  └─168.98° → 295850\n  └─172.65° → 296931\n  └─176.33° → 306528\n  └─180.00° → 372500"

  if visualtests
    @test_ref_plot "data/varioplane.png" plot(γ)
  end
end
