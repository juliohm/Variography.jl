# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVarioplane(data, var₁, var₂=var₁;
                        normal=spheredir(0,0),
                        nangs=50, ptol=0.5, dtol=0.5,
                        [parameters])

Given a `normal` direction, estimate the (cross-)variogram of variables
`var₁` and `var₂` along all directions in the corresponding plane of variation.

Optionally, specify the tolerance `ptol` for the plane partition, the tolerance
`dtol` for the direction partition, the number of angles `nangs` in the plane,
and forward the `parameters` to the underlying [`EmpiricalVariogram`](@ref).
"""
struct EmpiricalVarioplane{T,V}
  θs::Vector{T}
  γs::Vector{V}
end

function EmpiricalVarioplane(
  data::AbstractGeoTable,
  var₁,
  var₂=var₁;
  normal=spheredir(0, 0),
  nangs=50,
  ptol=0.5,
  dtol=0.5,
  kwargs...
)
  # sanity checks
  @assert nangs > 1 "nangs must be greater than one"

  # deterministic results
  rng = MersenneTwister(123)

  Dim = embeddim(domain(data))

  # basis for variogram plane
  if Dim == 2
    planes = [data]
    u, v = Vec(1.0, 0.0), Vec(0.0, 1.0)
  elseif Dim == 3
    planes = partition(rng, data, PlanePartition(normal, tol=ptol))
    u, v = householderbasis(normal)
  else
    @error "varioplane only supported in 2D or 3D"
  end

  # loop over half of the plane
  θs = range(0, stop=π, length=nangs)
  γs = map(θs) do θ
    dir = DirectionPartition(cos(θ) * u + sin(θ) * v, tol=dtol)

    γ(plane) = EmpiricalVariogram(partition(rng, plane, dir), var₁, var₂; kwargs...)
    foldxt(merge, Map(γ), collect(planes))
  end

  EmpiricalVarioplane(collect(θs), γs)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, ::EmpiricalVarioplane)
  print(io, "EmpiricalVarioplane")
end

function Base.show(io::IO, ::MIME"text/plain", γ::EmpiricalVarioplane)
  θs = [@sprintf "%.2f" rad2deg(θ) for θ in γ.θs]
  ns = [sum(values(g)[3]) for g in γ.γs]
  lines = ["  └─$(θ)° → $n" for (θ, n) in zip(θs, ns)]
  lines = length(lines) > 11 ? vcat(lines[1:5], ["  ⋮"], lines[(end - 4):end]) : lines
  println(io, γ)
  println(io, "  N° pairs")
  print(io, join(lines, "\n"))
end
