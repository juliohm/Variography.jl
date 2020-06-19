@testset "Empirical" begin
  # homogeneous field has zero variogram
  sdata = PointSetData(OrderedDict(:z=>ones(3)), Matrix(1.0I, 3, 3))
  γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2.)
  x, y, n = values(γ)
  @test x ≈ [1/2, 3/2]
  @test y[2] == 0.
  @test n == [0, 3]

  # basic test on number of lags
  sdata = PointSetData(OrderedDict(:z=>[1.,0.,1.]), [25. 50. 75.; 25. 75. 50.])
  γ = EmpiricalVariogram(sdata, :z, nlags=20, maxlag=1.)
  x, y, n = values(γ)
  @test length(x) == 20
  @test length(y) == 20
  @test length(n) == 20

  # empirical variogram on integer coordinates
  sdata = PointSetData(OrderedDict(:z=>ones(3)), Matrix(1I, 3, 3))
  γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2, algo=:full)
  x, y, n = values(γ)
  @test x ≈ [1/2, 3/2]
  @test y[2] == 0.
  @test n == [0, 3]

  # empirical variogram with only missing data
  X = rand(2,3)
  for z in [[NaN,NaN,NaN], [missing,missing,missing]]
    sdata = PointSetData(OrderedDict(:z=>z), X)
    γ = EmpiricalVariogram(sdata, :z, maxlag=1., nlags=5)
    x, y, n = values(γ)
    @test x == [.1, .3, .5, .7, .9]
    @test all(iszero.(n))
  end

  # accumulation algorithms give the same result
  sdata = PointSetData(OrderedDict(:z=>rand(1000)), rand(3,1000))
  γ₁ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algo=:full)
  γ₂ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algo=:ball)
  @test isequal(values(γ₁), values(γ₂))

  # print methods
  Random.seed!(123)
  d = RegularGridData(OrderedDict(:z=>rand(100,100)))
  γ = EmpiricalVariogram(d, :z)
  @test sprint(show, γ) == "EmpiricalVariogram"
  @test sprint(show, MIME"text/plain"(), γ) == "EmpiricalVariogram\n  abscissa: (0.35001785668734103, 13.650696410806301)\n  ordinate: (0.0, 0.083920131066808)\n  N° pairs: 2706158\n"

  if visualtests
    TI = training_image("WalkerLake")[1:20,1:20,1]
    d = RegularGridData(OrderedDict(:z=>TI))
    γ = EmpiricalVariogram(d, :z, maxlag=15.)
    @plottest plot(γ) joinpath(datadir,"empirical.png") !istravis
  end
end
