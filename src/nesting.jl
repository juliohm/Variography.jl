# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
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
struct NestedVariogram{CS,GS} <: Variogram
  cs::CS
  γs::GS

  function NestedVariogram{CS,GS}(cs, γs) where {CS,GS}
    @assert all(issymmetric.(cs)) "coefficients must be symmetric"
    new(cs, γs)
  end
end

NestedVariogram(cs, γs) = NestedVariogram{typeof(cs),typeof(γs)}(cs, γs)

(g::NestedVariogram)(h)                  = raw(sum(g.cs .* map(γ -> γ(h), g.γs)))
(g::NestedVariogram)(u::Point, v::Point) = raw(sum(g.cs .* map(γ -> γ(u, v), g.γs)))

sill(g::NestedVariogram)         = raw(sum(g.cs .* map(sill, g.γs)))
nugget(g::NestedVariogram)       = raw(sum(g.cs .* map(nugget, g.γs)))
Base.range(g::NestedVariogram)   = maximum(range(γ) for γ in g.γs)
isstationary(g::NestedVariogram) = all(isstationary(γ) for γ in g.γs)

# algebraic structure
*(c, γ::Variogram)                          = NestedVariogram((c,), (γ,))
*(c, γ::NestedVariogram)                    = NestedVariogram(map(x->c.*x, γ.cs), γ.γs)
+(γ₁::Variogram, γ₂::Variogram)             = NestedVariogram((1, 1), (γ₁, γ₂))
+(γ₁::NestedVariogram, γ₂::Variogram)       = NestedVariogram((γ₁.cs..., 1), (γ₁.γs..., γ₂))
+(γ₁::Variogram, γ₂::NestedVariogram)       = NestedVariogram((1, γ₂.cs...), (γ₁, γ₂.γs...))
+(γ₁::NestedVariogram, γ₂::NestedVariogram) = NestedVariogram((γ₁.cs..., γ₂.cs...), (γ₁.γs..., γ₂.γs...))

"""
    structures(γ)

Return the individual structures of a (possibly nested)
variogram as a tuple. The structures are the total nugget
`cₒ`, and the coefficients (or contributions) `cs` for the
remaining non-trivial structures `γs` after normalization
(i.e. sill=1, nugget=0).
"""
function structures(γ::Variogram)
  cₒ = nugget(γ)
  c  = sill(γ) - nugget(γ)
  T  = typeof(c)
  γ  = @set γ.sill = one(T)
  γ  = @set γ.nugget = zero(T)
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

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, g::NestedVariogram)
  N = length(g.cs)
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
