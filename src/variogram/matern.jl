# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MaternVariogram(range=r, sill=s, nugget=n, order=ν)
    MaternVariogram(ball; sill=s, nugget=n, order=ν)

A Matérn variogram with range `r`, sill `s` and nugget `n`.
The parameter `ν` is the order of the Bessel function.
Optionally, use a custom metric `ball`.
"""
struct MaternVariogram{V,O,B} <: Variogram
  sill::V
  nugget::V
  order::O
  ball::B
end

MaternVariogram(ball; sill=1.0, nugget=zero(typeof(sill)), order=1.0) = MaternVariogram(sill, nugget, order, ball)

MaternVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill)), order=1.0) =
  MaternVariogram(sill, nugget, order, MetricBall(range))

function (γ::MaternVariogram)(h::T) where {T}
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  ν = γ.order

  # shift lag by machine precision to
  # avoid explosion at the origin
  h′ = h + eps(T)

  δ = √(2ν) * 3(h′ / r)
  Β = besselk(ν, δ)
  Γ = gamma(ν)

  (s - n) * (1 - 2^(1 - ν) / Γ * δ^ν * Β) + (h > 0) * n
end

variotype(::MaternVariogram) = MaternVariogram

isstationary(::Type{<:MaternVariogram}) = true
