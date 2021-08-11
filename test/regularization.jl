@testset "Regularization" begin
  γ = GaussianVariogram()
  u = Point(0., 0.)
  v = Point(1., 0.)
  U = Quadrangle((-0.5,-0.5), (0.5,-0.5), (0.5,0.5), (-0.5,0.5))
  V = Quadrangle((0.5,-0.5), (1.5,-0.5), (1.5,0.5), (0.5,0.5))
  @test 0 < γ(U, v) < γ(u, v) < 1
  @test 0 < γ(u, V) < γ(u, v) < 1
  @test isapprox(γ(U, v), γ(u, V), atol=1e-1)
  @test 0 < γ(U, V) < γ(U, v) < 1
  @test 0 < γ(U, V) < γ(u, V) < 1
end
