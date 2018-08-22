@testset "Pairwise evaluation" begin
  Γ = pairwise(GaussianVariogram(), Matrix(1.0I, 3, 3))
  @test eltype(Γ) == Float64
  @test issymmetric(Γ)

  Γ_f = pairwise(GaussianVariogram(range=1f0, sill=1f0, nugget=0f0), Matrix(1.0f0I, 3, 3))
  @test eltype(Γ_f) == Float32
  @test issymmetric(Γ_f)

  grid = RegularGrid{Float64}(10, 10)
  Γ = pairwise(GaussianVariogram(), grid, 1:5)
  @test size(Γ) == (5, 5)
  @test issymmetric(Γ)
  Γ = pairwise(GaussianVariogram(), grid, 1:3, 7:10)
  @test size(Γ) == (3, 4)
  @test all(Γ .> 0)
end
