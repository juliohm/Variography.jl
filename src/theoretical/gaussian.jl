# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianVariogram(range=r, sill=s, nugget=n)
    GaussianVariogram(ball; sill=s, nugget=n)

A Gaussian variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct GaussianVariogram{T,B} <: Variogram
  sill::T
  nugget::T
  ball::B
end

GaussianVariogram(ball; sill=1.0, nugget=0.0) =
  GaussianVariogram(sill, nugget, ball)

GaussianVariogram(; range=1.0, sill=1.0, nugget=0.0) =
  GaussianVariogram(sill, nugget, MetricBall(range))

function (γ::GaussianVariogram)(h)
  # add small eps to nugget
  # for numerical stability
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget + typeof(s)(1e-6)
  (s - n) * (1 - exp(-3(h/r)^2)) + (h > 0) * n
end

isstationary(::Type{<:GaussianVariogram}) = true
