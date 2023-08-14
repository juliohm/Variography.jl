@testset "Sampling" begin
  γ = GaussianVariogram()
  seg = Segment((0.0, 0.0), (1.0, 1.0))
  ps = variosample(γ, seg) |> collect
  @test all(p -> Point(0.0, 0.0) ⪯ p ⪯ Point(1.0, 1.0), ps)
  @test length(ps) == 3

  γ = GaussianVariogram()
  quad = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
  ps = variosample(γ, quad) |> collect
  @test all(p -> Point(0.0, 0.0) ⪯ p ⪯ Point(1.0, 1.0), ps)
  @test length(ps) == 3 * 3

  γ = GaussianVariogram()
  hex = Hexahedron(
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0)
  )
  ps = variosample(γ, hex) |> collect
  @test all(p -> Point(0.0, 0.0, 0.0) ⪯ p ⪯ Point(1.0, 1.0, 1.0), ps)
  @test length(ps) == 3 * 3 * 3

  # regularization properties
  γ = GaussianVariogram()
  u = Point(0.0, 0.0)
  v = Point(1.0, 0.0)
  U = Quadrangle((-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5))
  V = Quadrangle((0.5, -0.5), (1.5, -0.5), (1.5, 0.5), (0.5, 0.5))
  @test 0 < γ(U, v) < γ(u, v) < 1
  @test 0 < γ(u, V) < γ(u, v) < 1
  @test isapprox(γ(U, v), γ(u, V), atol=1e-1)
  @test 0 < γ(U, V) < γ(U, v) < 1
  @test 0 < γ(U, V) < γ(u, V) < 1

  # deterministic samples in arbitrary geometries
  γ = GaussianVariogram()
  G = PolyArea((0.0, 0.0), (0.5, -1.5), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 1.5), (-0.5, 0.5))
  ps1 = variosample(γ, G)
  ps2 = variosample(γ, G)
  @test ps1 == ps2

  # samples with nugget effect model
  γ = NuggetEffect()
  h = Hexahedron(
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0)
  )
  ps = variosample(γ, h) |> collect
  @test length(ps) == 3 * 3 * 3
end
