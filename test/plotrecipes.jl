@testset "Plotting" begin
  if visualtests
    @testset "h-scatter" begin
      @plottest begin
        sdata = readgeotable(joinpath(datadir,"samples2D.tsv"), coordnames=(:x,:y))
        p0 = hscatter(sdata, :value, lag=0)
        p1 = hscatter(sdata, :value, lag=1)
        p2 = hscatter(sdata, :value, lag=2)
        p3 = hscatter(sdata, :value, lag=3)
        plot(p0, p1, p2, p3, layout=(2,2), size=(600,600))
      end joinpath(datadir,"hscatter.png") !istravis
    end
  end
end
