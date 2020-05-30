@testset "Varioplane" begin
  Z = readdlm(joinpath(datadir,"anisotropic.tsv"))
  sdata = RegularGridData{Float64}(OrderedDict(:Z => Z))
  vp = Varioplane(sdata, :Z, maxlag=50.)

  if visualtests
    # UNCOMMENT WHEN GR BACKEND ADDS SUPPORT FOR POLAR PLOTS
    # @plottest plot(vp) joinpath(datadir,"VarioPlane.png") !travis
  end
end
