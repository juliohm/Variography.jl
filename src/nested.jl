# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NestedVariogram(cs, γs)

A nested variogram model `γ = c₁γ₁ + c₂γ₂ + ⋯ + cₙγₙ` with
coefficients `cs = (c₁, c₂, ..., cₙ)` and variogram models
`γs = (γ₁, γ₂, ..., γₙ)`.
"""
struct NestedVariogram{N} <: Variogram{Number,Metric}
  cs::NTuple{N,Any}
  γs::NTuple{N,Variogram}
end
(g::NestedVariogram)(h) = sum(c*γ(h) for (c, γ) in zip(g.cs, g.γs))
(g::NestedVariogram)(x, y) = sum(c*γ(x,y) for (c, γ) in zip(g.cs, g.γs))
Base.range(g::NestedVariogram) = maximum(range(γ) for γ in g.γs)
sill(g::NestedVariogram) = sum(c*sill(γ) for (c, γ) in zip(g.cs, g.γs))
nugget(g::NestedVariogram) = sum(c*nugget(γ) for (c, γ) in zip(g.cs, g.γs))
isstationary(g::NestedVariogram) = all(isstationary(γ) for γ in g.γs)

# algebraic structure
*(c, γ::Variogram) = NestedVariogram((c,), (γ,))
*(c, γ::NestedVariogram{N}) where {N} = NestedVariogram(ntuple(i->c.*γ.cs[i], N), γ.γs)
+(γ₁::Variogram, γ₂::Variogram) = NestedVariogram((I, I), (γ₁, γ₂))
+(γ₁::NestedVariogram, γ₂::Variogram) = NestedVariogram((γ₁.cs..., I), (γ₁.γs..., γ₂))
+(γ₁::Variogram, γ₂::NestedVariogram) = NestedVariogram((I, γ₂.cs...), (γ₁, γ₂.γs...))
+(γ₁::NestedVariogram, γ₂::NestedVariogram) = NestedVariogram((γ₁.cs..., γ₂.cs...), (γ₁.γs..., γ₂.γs...))

# ------------
# IO methods
# ------------
function Base.show(io::IO, ::NestedVariogram{N}) where N
  print(io, "NestedVariogram{$N}")
end

function Base.show(io::IO, ::MIME"text/plain", g::NestedVariogram)
  println(io, g)
  println(io, "  models")
  coeffs = [c isa UniformScaling ? 1.0 : c for c in g.cs]
  models = [nameof(typeof(γ)) for γ in g.γs]
  params = ["range=$(range(γ)), sill=$(sill(γ)), nugget=$(nugget(γ))" for γ in g.γs]
  lines  = ["    └─$γ($p)" for (γ, p) in zip(models, params)]
  println(io, join(lines, "\n"))
  println(io, "  coefficients")
  lines  = ["    └─$c" for c in coeffs]
  print(io, join(lines, "\n"))
end
