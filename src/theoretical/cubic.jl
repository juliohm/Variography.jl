# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CubicVariogram(range=r, sill=s, nugget=n)
    CubicVariogram(ball; sill=s, nugget=n)

A cubic variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct CubicVariogram{V,B} <: Variogram
  sill::V
  nugget::V
  ball::B
end

CubicVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = CubicVariogram(sill, nugget, ball)

CubicVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) = CubicVariogram(sill, nugget, MetricBall(range))

function (γ::CubicVariogram)(h::T) where {T}
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget

  # constants
  c1 = T(35) / T(4)
  c2 = T(7) / T(2)
  c3 = T(3) / T(4)
  s1 = 7 * (h / r)^2 - c1 * (h / r)^3 + c2 * (h / r)^5 - c3 * (h / r)^7
  s2 = T(1)

  (h < r) * (s - n) * s1 + (h ≥ r) * (s - n) * s2 + (h > 0) * n
end

variotype(::CubicVariogram) = CubicVariogram

isstationary(::Type{<:CubicVariogram}) = true
