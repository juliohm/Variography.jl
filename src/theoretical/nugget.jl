# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NuggetEffect(nugget=n)

A pure nugget effect variogram with nugget `n`.
"""
struct NuggetEffect{T} <: Variogram
  nugget::T
end

NuggetEffect(; nugget=1.0) = NuggetEffect(nugget)

(γ::NuggetEffect)(h) = (h > 0) * γ.nugget

(γ::NuggetEffect)(u::Point, v::Point) =
  γ(evaluate(Euclidean(), coordinates(u), coordinates(v)))

Base.range(::NuggetEffect{T}) where {T} = zero(T)

sill(γ::NuggetEffect) = γ.nugget

isstationary(::Type{<:NuggetEffect}) = true