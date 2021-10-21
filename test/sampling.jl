@testset "Sampling" begin
  γ = GaussianVariogram()
  seg = Segment((0.,0.),(1.,1.))
  ps = Variography._sample(γ, seg)
  @test all(p -> Point(0.,0.) ⪯ p ⪯ Point(1.,1.), ps)
  @test length(ps) == 3

  γ = GaussianVariogram()
  quad = Quadrangle((0.,0.), (1.,0.), (1.,1.), (0.,1.))
  ps = Variography._sample(γ, quad)
  @test all(p -> Point(0.,0.) ⪯ p ⪯ Point(1.,1.), ps)
  @test length(ps) == 3*3

  γ = GaussianVariogram()
  hex = Hexahedron((0.,0.,0.), (1.,0.,0.), (1.,1.,0.), (0.,1.,0.),
                   (0.,0.,1.), (1.,0.,1.), (1.,1.,1.), (0.,1.,1.))
  ps = Variography._sample(γ, hex)
  @test all(p -> Point(0.,0.,0.) ⪯ p ⪯ Point(1.,1.,1.), ps)
  @test length(ps) == 3*3*3

  # regularization properties
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

  # deterministic samples in arbitrary geometries
  γ = GaussianVariogram()
  G = PolyArea((0.,0.), (0.5,-1.5), (1.,0.), (1.5,0.5),
               (1.,1.), (0.5,1.5), (-0.5,0.5), (0.,0.))
  ps1 = Variography._sample(γ, G)
  ps2 = Variography._sample(γ, G)
  @test ps1 == ps2
end
