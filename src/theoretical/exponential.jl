# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialVariogram(range=r, sill=s, nugget=n)
    ExponentialVariogram(ball; sill=s, nugget=n)

An exponential variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct ExponentialVariogram{T,B} <: Variogram
  sill::T
  nugget::T
  ball::B
end

ExponentialVariogram(ball; sill=1.0, nugget=0.0) =
  ExponentialVariogram(sill, nugget, ball)

ExponentialVariogram(; range=1.0, sill=1.0, nugget=0.0) =
  ExponentialVariogram(sill, nugget, MetricBall(range))

function (γ::ExponentialVariogram)(h)
  r = boundaryvalue(γ.ball)
  s = γ.sill
  n = γ.nugget
  (s - n) * (1 - exp(-3(h/r))) + (h > 0) * n
end

isstationary(::Type{<:ExponentialVariogram}) = true
