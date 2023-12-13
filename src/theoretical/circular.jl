"""
    CircularVariogram(; range=r, sill=s, nugget=n)
    CircularVariogram(ball; sill=s, nugget=n)

A circular variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct CircularVariogram{V,B} <: Variogram
  sill::V
  nugget::V
  ball::B
end

CircularVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = CircularVariogram(sill, nugget, ball)

CircularVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) = CircularVariogram(sill, nugget, MetricBall(range))

function (γ::CircularVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  v = h ≤ r ? 1 - (2 / π) * acos(h / r) + (2h / (π * r)) * sqrt(1 - (h^2 / r^2)) : one(h)
  (s - n) * v + (h > zero(h)) * n
end

variotype(::CircularVariogram) = CircularVariogram

isstationary(::Type{<:CircularVariogram}) = true
