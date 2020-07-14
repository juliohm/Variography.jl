@testset "Plotting" begin
  if visualtests
    @testset "h-scatter" begin
      @plottest begin
        sdata = readgeotable(joinpath(datadir,"samples2D.tsv"), coordnames=(:x,:y))
        hscatter(sdata, :value, lags=[0.,1.,2.,3.], layout=(2,2), size=(600,600))
      end joinpath(datadir,"hscatter.png") !istravis
    end
  end
end
