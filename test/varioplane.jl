@testset "Varioplane" begin
  Random.seed!(123)
  Z = readdlm(joinpath(datadir,"anisotropic.tsv"))
  sdata = RegularGridData{Float64}(OrderedDict(:Z => Z))
  vp = EmpiricalVarioplane(sdata, :Z, maxlag=50.)

  if visualtests
    @plottest plot(vp) joinpath(datadir,"varioplane.png") !istravis
  end
end
