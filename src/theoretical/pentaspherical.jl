# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PentasphericalVariogram(range=r, sill=s, nugget=n)
    PentasphericalVariogram(ball; sill=s, nugget=n)

A pentaspherical variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct PentasphericalVariogram{T,B} <: Variogram
  sill::T
  nugget::T
  ball::B
end

PentasphericalVariogram(ball; sill=1.0, nugget=0.0) =
  PentasphericalVariogram(sill, nugget, ball)

PentasphericalVariogram(; range=1.0, sill=1.0, nugget=0.0) =
  PentasphericalVariogram(sill, nugget, MetricBall(range))

function (γ::PentasphericalVariogram)(h::T) where {T}
  r = boundaryvalue(γ.ball)
  s = γ.sill
  n = γ.nugget

  # constants
  c1 = T(15) / T(8)
  c2 = T(5)  / T(4)
  c3 = T(3)  / T(8)
  s1 = c1*(h/r) - c2*(h/r)^3 + c3*(h/r)^5
  s2 = T(1)

  (h < r) * (s - n) * s1 +
  (h ≥ r) * (s - n) * s2 +
  (h > 0) * n
end

isstationary(::Type{<:PentasphericalVariogram}) = true
