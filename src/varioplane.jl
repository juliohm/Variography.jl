# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVarioplane(sdata, var₁, var₂=var₁;
                        normal=spheredir(0,0), ptol=1e-6, dtol=0.5,
                        nangs=50, nlags=20, maxlag=nothing)

Given a `normal` direction, estimate the (cross-)variogram of variables
`var₁` and `var₂` along all directions in the corresponding plane of variation.

Optionally, specify the tolerance `ptol` for the [`PlanePartitioner`](@ref)
the tolerance `dtol` for the [`DirectionPartitioner`](@ref), the number of
angles `nangs` in the plane, the number of lags `nlags` along each of these
angles, and the maximum lag `maxlag` to consider.
"""
struct EmpiricalVarioplane{T,V}
  θs::Vector{T}
  γs::Vector{V}
end

function EmpiricalVarioplane(sdata, var₁::Symbol, var₂::Symbol=var₁;
                             normal=spheredir(0,0), ptol=1e-6, dtol=0.5,
                             nangs=50, nlags=20, maxlag=nothing)
  # sanity checks
  @assert nangs > 1 "nangs must be greater than one"

  # basis for variogram plane
  if ndims(sdata) == 2
    planes = [sdata]
    u, v = SVector(1.,0.), SVector(0.,1.)
  elseif ndims(sdata) == 3
    planes = partition(sdata, PlanePartitioner(normal, tol=ptol))
    u, v = planebasis(normal)
  else
    @error "varioplane only supported in 2D or 3D"
  end

  # loop over half of the plane
  θs = range(0, stop=π, length=nangs)
  γs = map(θs) do θ
    dir = DirectionPartitioner(cos(θ)*u + sin(θ)*v, tol=dtol)

    # compute directional variogram across planes
    function γ(plane)
      p = partition(plane, dir)
      EmpiricalVariogram(p, var₁, var₂; nlags=nlags, maxlag=maxlag)
    end
    reduce(merge, Map(γ), collect(planes))
  end

  EmpiricalVarioplane(collect(θs), γs)
end

# ------------
# IO methods
# ------------
function Base.show(io::IO, γ::EmpiricalVarioplane)
  print(io, "EmpiricalVarioplane")
end

function Base.show(io::IO, ::MIME"text/plain", γ::EmpiricalVarioplane)
  θs = [@sprintf "%.2f" rad2deg(θ) for θ in γ.θs]
  ns = [sum(values(g)[3]) for g in γ.γs]
  lines = ["  └─$(θ)° → $n" for (θ, n) in zip(θs, ns)]
  lines = length(lines) > 11 ? vcat(lines[1:5],["  ⋮"],lines[end-4:end]) : lines
  println(io, γ)
  println(io, "  N° pairs")
  print(io, join(lines, "\n"))
end
