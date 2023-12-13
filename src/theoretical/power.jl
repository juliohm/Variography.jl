# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PowerVariogram(scaling=s, exponent=a, nugget=n)

A power variogram with scaling `s`, exponent `a` and nugget `n`.
"""
struct PowerVariogram{V,E} <: Variogram
  scaling::V
  nugget::V
  exponent::E
end

PowerVariogram(; scaling=1.0, nugget=zero(typeof(scaling)), exponent=1.0) = PowerVariogram(scaling, nugget, exponent)

function (γ::PowerVariogram)(h)
  s = γ.scaling
  a = γ.exponent
  n = γ.nugget
  s * h^a + (h > 0) * n
end

function (γ::PowerVariogram)(u::Point, v::Point)
  d = Euclidean()
  x = coordinates(u)
  y = coordinates(v)
  h = evaluate(d, x, y)
  γ(h)
end

variotype(::PowerVariogram) = PowerVariogram

isstationary(::Type{<:PowerVariogram}) = false

isisotropic(::PowerVariogram) = true
