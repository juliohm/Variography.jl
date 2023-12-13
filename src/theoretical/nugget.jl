# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NuggetEffect(nugget=n)

A pure nugget effect variogram with nugget `n`.
"""
struct NuggetEffect{V} <: Variogram
  nugget::V
end

NuggetEffect(; nugget=1.0) = NuggetEffect(nugget)

(γ::NuggetEffect)(h) = (h > 0) * γ.nugget

function (γ::NuggetEffect)(u::Point, v::Point)
  d = Euclidean()
  x = coordinates(u)
  y = coordinates(v)
  h = evaluate(d, x, y)
  γ(h)
end

Base.range(::NuggetEffect{T}) where {T} = zero(T)

sill(γ::NuggetEffect) = γ.nugget

variotype(::NuggetEffect) = NuggetEffect

isstationary(::Type{<:NuggetEffect}) = true

isisotropic(::NuggetEffect) = true
