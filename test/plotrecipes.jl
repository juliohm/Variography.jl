@testset "Plotting" begin
  if ismaintainer || istravis
    @testset "h-scatter" begin
      function plot_hscatter(fname)
        hscatter(geodf2D, :value, lags=[0.,1.,2.,3.], size=(1000,300))
        png(fname)
      end
      refimg = joinpath(datadir,"HScatter.png")
      @test test_images(VisualTest(plot_hscatter, refimg), popup=!istravis, tol=0.15) |> success
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
