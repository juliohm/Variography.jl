# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Varioplane(sdata, var₁, var₂=var₁,
               theta=0, phi=90, ptol=1e-6, dtol=0.5,
               nangs=50, nlags=20, maxlag=nothing)

Given a normal direction in spherical coordinates `theta` and `phi` in degrees,
estimate the (cross-)variogram of variables `var₁` and `var₂` along all directions
in the corresponding plane of variation.

Optionally, specify the tolerance `ptol` for the [`PlanePartitioner`](@ref)
the tolerance `dtol` for the [`DirectionPartitioner`](@ref), the number of
angles `nangs` in the plane, the number of lags `nlags` along each of these
angles, and the maximum lag `maxlag` to consider.
"""
struct Varioplane{T,V<:EmpiricalVariogram}
  θs::Vector{T}
  γs::Vector{V}
end

function Varioplane(sdata, var₁::Symbol, var₂::Symbol=var₁;
                    theta=0, phi=90, ptol=1e-6, dtol=0.5,
                    nangs=50, nlags=20, maxlag=nothing)
  # sanity checks
  @assert 0 ≤ theta ≤ 360 "theta must lie in [0,360]"
  @assert -90 ≤ phi ≤ 90 "phi must lie in [-90,90]"
  @assert nangs > 1 "nangs must be greater than one"

  # basis for variogram plane
  N = ndims(sdata)
  if N == 2
    planes = [sdata]
    u, v = (1.,0.), (0.,1.)
  elseif N == 3
    θ = deg2rad(theta); ϕ = deg2rad(phi)
    n = (cos(ϕ)cos(θ), cos(ϕ)sin(θ), sin(ϕ))

    planes = partition(sdata, PlanePartitioner(n, tol=ptol))
    u, v = planebasis(n)
  else
    @error "varioplane only supported in 2D or 3D"
  end

  # loop over half of the plane
  θs = range(0, stop=π, length=nangs)
  γs = map(θs) do θ
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

    γ
  end

  Varioplane(collect(θs), γs)
end
