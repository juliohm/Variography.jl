@testset "Covariance" begin
  x, y = rand(Point3), rand(Point3)
  for (CovType, VarioType) in [
    (CircularCovariance, CircularVariogram),
    (CubicCovariance, CubicVariogram),
    (ExponentialCovariance, ExponentialVariogram),
    (GaussianCovariance, GaussianVariogram),
    (MaternCovariance, MaternVariogram),
    (NuggetCovariance, NuggetEffect),
    (PentasphericalCovariance, PentasphericalVariogram),
    (SineHoleCovariance, SineHoleVariogram),
    (SphericalCovariance, SphericalVariogram)
  ]
    γ = VarioType()
    cov = CovType(γ)
    @test cov(x, y) == sill(γ) - γ(x, y)
  end

  # shows
  cov = CircularCovariance()
  @test sprint(show, cov) == "CircularCovariance(sill: 1.0, nugget: 0.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  CircularCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0
  └─ distance: Euclidean"""

  cov = CubicCovariance()
  @test sprint(show, cov) == "CubicCovariance(sill: 1.0, nugget: 0.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  CubicCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0
  └─ distance: Euclidean"""

  cov = ExponentialCovariance()
  @test sprint(show, cov) == "ExponentialCovariance(sill: 1.0, nugget: 0.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  ExponentialCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0
  └─ distance: Euclidean"""

  cov = GaussianCovariance()
  @test sprint(show, cov) == "GaussianCovariance(sill: 1.0, nugget: 0.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  GaussianCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0
  └─ distance: Euclidean"""

  cov = MaternCovariance()
  @test sprint(show, cov) == "MaternCovariance(sill: 1.0, nugget: 0.0, order: 1.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  MaternCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ order: 1.0
  ├─ range: 1.0
  └─ distance: Euclidean"""

  cov = NuggetCovariance()
  @test sprint(show, cov) == "NuggetCovariance(nugget: 1.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  NuggetCovariance
  └─ nugget: 1.0"""

  cov = PentasphericalCovariance()
  @test sprint(show, cov) == "PentasphericalCovariance(sill: 1.0, nugget: 0.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  PentasphericalCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0
  └─ distance: Euclidean"""

  cov = SineHoleCovariance()
  @test sprint(show, cov) == "SineHoleCovariance(sill: 1.0, nugget: 0.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  SineHoleCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0
  └─ distance: Euclidean"""

  cov = SphericalCovariance()
  @test sprint(show, cov) == "SphericalCovariance(sill: 1.0, nugget: 0.0, range: 1.0, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  SphericalCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0
  └─ distance: Euclidean"""
end
