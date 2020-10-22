# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper function to extract raw data
# from uniform scaling objects
raw(a::UniformScaling) = a.λ
raw(a) = a

"""
    NestedVariogram(cs, γs)

A nested variogram model `γ = c₁γ₁ + c₂γ₂ + ⋯ + cₙγₙ` with
coefficients `cs = (c₁, c₂, ..., cₙ)` and variogram models
`γs = (γ₁, γ₂, ..., γₙ)`.
"""
struct NestedVariogram{N} <: Variogram{Number,Metric}
  cs::NTuple{N,Any}
  γs::NTuple{N,Variogram}

  function NestedVariogram{N}(cs, γs) where {N}
    @assert all(issymmetric.(cs)) "coefficients must be symmetric"
    new(cs, γs)
  end
end

NestedVariogram(cs, γs) = NestedVariogram{length(cs)}(cs, γs)

# variogram interface
(g::NestedVariogram)(h)          = raw(sum(c*γ(h) for (c, γ) in zip(g.cs, g.γs)))
(g::NestedVariogram)(x, y)       = raw(sum(c*γ(x, y) for (c, γ) in zip(g.cs, g.γs)))
sill(g::NestedVariogram)         = raw(sum(c*sill(γ) for (c, γ) in zip(g.cs, g.γs)))
nugget(g::NestedVariogram)       = raw(sum(c*nugget(γ) for (c, γ) in zip(g.cs, g.γs)))
Base.range(g::NestedVariogram)   = maximum(range(γ) for γ in g.γs)
isstationary(g::NestedVariogram) = all(isstationary(γ) for γ in g.γs)

# algebraic structure
*(c, γ::Variogram)                          = NestedVariogram((c,), (γ,))
*(c, γ::NestedVariogram)                    = NestedVariogram(map(x->c.*x, γ.cs), γ.γs)
+(γ₁::Variogram, γ₂::Variogram)             = NestedVariogram((I, I), (γ₁, γ₂))
+(γ₁::NestedVariogram, γ₂::Variogram)       = NestedVariogram((γ₁.cs..., I), (γ₁.γs..., γ₂))
+(γ₁::Variogram, γ₂::NestedVariogram)       = NestedVariogram((I, γ₂.cs...), (γ₁, γ₂.γs...))
+(γ₁::NestedVariogram, γ₂::NestedVariogram) = NestedVariogram((γ₁.cs..., γ₂.cs...), (γ₁.γs..., γ₂.γs...))

"""
    structures(γ)

Return the individual structures of a (possibly nested)
variogram. The structures are the total nugget `cₒ`, and
the coefficients (or contributions) `cs` for the remaining
non-trivial structures `γs`.
"""
function structures(γ::Variogram)
  cₒ = nugget(γ)
  c  = sill(γ) - nugget(γ)
  T  = typeof(c)
  γ = @set γ.sill = one(T)
  γ = @set γ.nugget = zero(T)
  cₒ, (c,), (γ,)
end

function structures(γ::NestedVariogram)
  ks, gs = γ.cs, γ.γs

  # total nugget and contributions
  cₒ = raw(sum(@. ks * nugget(gs)))
  cs = @. raw(ks * (sill(gs) - nugget(gs)))

  # discard nugget effect terms
  inds = findall(g->!(g isa NuggetEffect), gs)
  cs, γs = cs[inds], gs[inds]

  # adjust sill and nugget
  γs = map(γs) do γ
    T = typeof(sill(γ))
    γ = @set γ.sill = one(T)
    γ = @set γ.nugget = zero(T)
  end

  cₒ, cs, γs
end

# ------------
# IO methods
# ------------
function Base.show(io::IO, ::NestedVariogram{N}) where N
  print(io, "NestedVariogram{$N}")
end

function Base.show(io::IO, ::MIME"text/plain", g::NestedVariogram)
  coeffs = 1 .* raw.(g.cs)
  models = [nameof(typeof(γ)) for γ in g.γs]
  params = ["range=$(range(γ)), sill=$(sill(γ)), nugget=$(nugget(γ))" for γ in g.γs]
  println(io, g)
  println(io, "  structures")
  lines  = ["    └─$γ($p)" for (γ, p) in zip(models, params)]
  println(io, join(lines, "\n"))
  println(io, "  coefficients")
  lines  = ["    └─$c" for c in coeffs]
  print(io, join(lines, "\n"))
end
