# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

@userplot VarPlane

@recipe function f(vp::VarPlane;
                   theta=0, phi=90,
                   ptol=1e-6, dtol=0.5,
                   nangs=50, nlags=20, maxlag=nothing,
                   showrange=true, varmodel=GaussianVariogram)
  # sanity checks
  @assert 0 ≤ theta ≤ 360 "theta must lie in [0,360]"
  @assert -90 ≤ phi ≤ 90 "phi must lie in [-90,90]"

  # get inputs
  spatialdata = vp.args[1]
  var₁ = vp.args[2]
  var₂ = length(vp.args) == 3 ? vp.args[3] : var₁
  N = ndims(spatialdata)

  # basis for variogram plane
  if N == 2
    planes = [spatialdata]
    u, v = (1.,0.), (0.,1.)
  elseif N == 3
    θ = deg2rad(theta); ϕ = deg2rad(phi)
    n = (cos(ϕ)cos(θ), cos(ϕ)sin(θ), sin(ϕ))

    planes = partition(spatialdata, PlanePartitioner(n, tol=ptol))
    u, v = plane_basis(n)
  else
    @error "variogram plane only supported in 2D or 3D"
  end

  # polar coordinates and variogram values
  rs = Vector{Float64}(undef, nlags)
  θs = range(0, stop=π, length=nangs)
  zs = Matrix{Float64}(undef, nlags, nangs)

  # ranges for each angle
  ls = Vector{Float64}(undef, nangs)

  # loop over half of the plane
  for (j, θ) in Iterators.enumerate(θs)
    dir = ntuple(i -> cos(θ)*u[i] + sin(θ)*v[i], N)

    dpart = DirectionPartitioner(dir, tol=dtol)

    # compute directional variogram across planes
    plane, _ = iterate(planes)
    γ = EmpiricalVariogram(partition(plane, dpart), var₁, var₂;
                           nlags=nlags, maxlag=maxlag)
    for plane in Iterators.drop(planes, 1)
      γplane = EmpiricalVariogram(partition(plane, dpart), var₁, var₂;
                                  nlags=nlags, maxlag=maxlag)
      merge!(γ, γplane)
    end

    rs, ys, _ = values(γ)

    # clean NaN values (i.e. empty bins)
    isnan(ys[1]) && (ys[1] = 0)
    for i in 2:length(ys)
      isnan(ys[i]) && (ys[i] = ys[i-1])
    end

    # save variogram for this angle
    zs[:,j] = ys

    if showrange
      γtheo = fit(varmodel, γ)
      ls[j] = range(γtheo)
    end
  end

  # exploit symmetry
  θs = range(0, stop=2π, length=2nangs)
  zs = hcat(zs, zs)

  # show variables in plot title
  ptitle = var₁ == var₂ ? "$var₁" : "$var₁-$var₂"

  projection --> :polar
  ylim --> (0, Inf)
  title --> ptitle

  # plot variogram plane
  @series begin
    seriestype --> :heatmap
    color --> :curl

    θs, rs, zs
  end

  # plot ranges
  if showrange
    @series begin
      seriestype --> :path
      linecolor --> :black
      label --> "range"

      θs, ls
    end
  end
end

function plane_basis(normal::NTuple{3,T}) where {T<:Real}
  # normalize input
  n = normal ./ sqrt(sum(normal[i]^2 for i in 1:3))

  # find last non-zero component
  idx = -1
  for (i, c) in enumerate(reverse(n))
    if c != 0
      idx = length(n) - i + 1
      break
    end
  end

  @assert idx > 0 "invalid normal vector"

  # first basis vector (perturb and subtract projection)
  u = ntuple(i -> i == idx%3 + 1 ? n[i] + one(T) : n[i], 3)
  l = sum(u[i]*n[i] for i in 1:3)
  u = ntuple(i -> u[i] - l*n[i], 3)

  # second basis vector (cross product)
  nx, ny, nz = n
  ux, uy, uz = u
  v = (ny*uz - nz*uz, nz*ux - nx*uz, nx*uy - ny*ux)

  # normalize output
  u = u ./ sqrt(sum(u[i]^2 for i in 1:3))
  v = v ./ sqrt(sum(v[i]^2 for i in 1:3))

  u, v
end
