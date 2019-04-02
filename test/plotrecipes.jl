@testset "Plotting" begin
  if visualtests
    @testset "h-scatter" begin
      @plottest begin
        hscatter(geodf2D, :value, lags=[0.,1.,2.,3.], size=(1000,300))
      end joinpath(datadir,"HScatter.png") !istravis
    end

    # UNCOMMENT WHEN GR ADDS SUPPORT TO POLAR PLOTS
    # @testset "varplane" begin
      # imgdata = readdlm(joinpath(datadir,"anisotropic.tsv"))
      # geodata = RegularGridData{Float64}(Dict(:z => imgdata))

      # function plot_varplane(fname)
        # varplane(geodata, :z, maxlag=50., size=(500,500))
        # png(fname)
      # end
      # refimg = joinpath(datadir,"VarPlane.png")
      # @test test_images(VisualTest(plot_varplane, refimg), popup=!istravis) |> success
    # end
  end
end
