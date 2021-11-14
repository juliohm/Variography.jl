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
struct MaternVariogram{T,B} <: Variogram
  sill::T
  nugget::T
  order::T
  ball::B
end

MaternVariogram(ball; sill=1.0, nugget=0.0, order=1.0) =
  MaternVariogram(sill, nugget, order, ball)

MaternVariogram(; range=1.0, sill=1.0, nugget=0.0, order=1.0) =
  MaternVariogram(sill, nugget, order, MetricBall(range))

function (γ::MaternVariogram)(h::T) where {T}
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  ν = γ.order

  # shift lag by machine precision to
  # avoid explosion at the origin
  h′ = sqrt(2ν)*(h + eps(T)) / r
  Β  = besselk(ν, h′)
  Γ  = gamma(ν)

  (s - n) * (1 - 2^(1 - ν) / Γ * h′^ν * Β) + (h′ > 0) * n
end

isstationary(::Type{<:MaternVariogram}) = true
