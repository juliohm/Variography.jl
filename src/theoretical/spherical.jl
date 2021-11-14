# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SphericalVariogram(range=r, sill=s, nugget=n)
    SphericalVariogram(ball; sill=s, nugget=n)

A spherical variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct SphericalVariogram{T,B} <: Variogram
  sill::T
  nugget::T
  ball::B
end

SphericalVariogram(ball; sill=1.0, nugget=0.0) =
  SphericalVariogram(sill, nugget, ball)

SphericalVariogram(; range=1.0, sill=1.0, nugget=0.0) =
  SphericalVariogram(sill, nugget, MetricBall(range))

function (γ::SphericalVariogram)(h::T) where {T}
  r = boundaryvalue(γ.ball)
  s = γ.sill
  n = γ.nugget

  # constants
  c1 = T(3) / T(2)
  c2 = T(1) / T(2)
  s1 = c1*(h/r) - c2*(h/r)^3
  s2 = T(1)

  (h < r) * (s - n) * s1 +
  (h ≥ r) * (s - n) * s2 +
  (h > 0) * n
end

isstationary(::Type{<:SphericalVariogram}) = true
