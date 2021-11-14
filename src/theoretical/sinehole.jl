# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SineHoleVariogram(range=r, sill=s, nugget=n)
    SineHoleVariogram(ball; sill=s, nugget=n)

A sine hole variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct SineHoleVariogram{T,B} <: Variogram
  sill::T
  nugget::T
  ball::B
end

SineHoleVariogram(ball; sill=1.0, nugget=0.0) =
  SineHoleVariogram(sill, nugget, ball)

SineHoleVariogram(; range=1.0, sill=1.0, nugget=0.0) =
  SineHoleVariogram(sill, nugget, MetricBall(range))

function (γ::SineHoleVariogram)(h::T) where {T}
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget

  # shift lag by machine precision to
  # avoid explosion at the origin
  h′ = h + eps(T)
  c  = T(π)

  (s - n) * (1 - sin(c*h′/r)/(c*h′/r)) + (h′ > 0) * n
end

isstationary(::Type{<:SineHoleVariogram}) = true
