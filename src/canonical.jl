# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NestedVariogram{T,D}

A theoretical variogram model with nugget effect and nested strucutures
with parameters of type `T` and distance of type `D`.
"""
@with_kw struct NestedVariogram{T,D} <: Variogram{T,D}
  nugget::T = 0.0
  scales::Array{T} = T[]
  structures::Array{Variogram} = Variogram[]
  distance::D = Euclidean()
end

"""
    range(γ)

Return the range of the variogram `γ` when defined.
"""
Base.range(γ::NestedVariogram) = maximum([range(st) for st in γ.structures])

"""
    distance(γ)

Return the distance of the variogram `γ`.
"""
distance(γ::NestedVariogram) = γ.distance

"""
    γ(h)

Evaluate the variogram at distance `h`.
"""
(γ::NestedVariogram)(h) = (h > 0) * γ.nugget + sum([st(h) for st in γ.structures] .* γ.scales)

sill(γ::NestedVariogram) = γ.nugget + sum(γ.scales)
isstationary(::Type{<:NestedVariogram}) = prod([isstationary(st) for st in γ.structures])

"""
    γ₁ + γ₂

Fuse two any Variogram `γ₁` and `γ₂` by fusing transformed ones to canonical.
"""
+(γ₁::Variogram, γ₂::Variogram) = canonical(γ₁) + canonical(γ₂)

"""
    γ₁ + γ₂

Fuse two NestedVariogram `γ₁` and `γ₂`.
"""
function +(γ₁::NestedVariogram, γ₂::NestedVariogram)
  nugget_fused = γ₁.nugget + γ₂.nugget
  scales_fused = vcat(γ₁.scales,γ₂.scales)
  structures_fused = vcat(γ₁.structures, γ₂.structures)
  γ_nested = NestedVariogram(nugget=nugget_fused, scales=scales_fused, structures=structures_fused)
  # TODO take care of distance
  return γ_nested
end

"""
    c * γ

Scale any Variogram `γ` by a constant `c` by scaling the transformed one to canonical.
"""
*(c::Number, γ::Variogram) = c*canonical(γ)

"""
    c * γ

Scale a NestedVariogram `γ` by a constant `c` by scaling inner nugget and scales.
"""
*(c::Number, γ::NestedVariogram) = NestedVariogram(
      nugget=γ.nugget*c,
      scales=γ.scales*c,
      structures=γ.structures,
      distance=γ.distance
  )

"""
    canonical(γ)

Return for any variogram `γ` its canonical form as a nested variogram.
"""
function canonical(γ::Variogram)
  if !isa(γ, NestedVariogram)
    nugget_comp, scales, structures = extract_variogram(γ)
    γ_nested = NestedVariogram(nugget=nugget_comp, scales=scales, structures=structures,distance=Euclidean())

    # TODO check if all distances are the same to take advantage of

    γ_nested
  else
    γ
  end
end

#-------------------------------
# CANONICAL STRUCTURE CONVERSION
#-------------------------------
canonical_structure(γ::SphericalVariogram{T,D}) where {T,D} = 
    SphericalVariogram(sill=one(T),nugget=zero(T), distance=γ.distance, range=γ.range)

canonical_structure(γ::ExponentialVariogram{T,D}) where {T,D} = 
    ExponentialVariogram(sill=one(T),nugget=zero(T), distance=γ.distance, range=γ.range)

"""
    extract_variogram(γ)

Return for any variogram `γ` its nugget, scales and canonical structures.
"""
function extract_variogram(γ::Variogram{T,D}) where {T,D}
  scales = T[]
  structures = Variogram[]
  if isa(γ, SumVariogram)
    nugget1, scales1, structs1 = extract_variogram(γ.γ1)
    nugget2, scales2, structs2 = extract_variogram(γ.γ2)

    nugget_comp = nugget1 + nugget2
    append!(scales, scales1)
    append!(scales, scales2)
    append!(structures, structs1)
    append!(structures, structs2)
  elseif isa(γ, ScaledVariogram)
    sill_scaled = γ.c
    nugget_comp, scales_inner, structs_inner = extract_variogram(γ.γ)
    nugget_comp = nugget_comp * sill_scaled
    append!(scales, scales_inner*sill_scaled)
    append!(structures, structs_inner)
  else
    sill_comp = sill(γ)
    nugget_comp = nugget(γ)
    if sill_comp > 0
      push!(scales, sill_comp - nugget_comp)
      # redefine the structure to its canonical form (sill=1.0 and nugget=0.0)
      γ_structure = canonical_structure(γ)
      push!(structures, γ_structure)
    end
  end

  return nugget_comp, scales, structures
end
