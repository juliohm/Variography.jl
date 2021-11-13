# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Variogram{T,D}

A theoretical variogram model (e.g. Gaussian variogram)
with parameters of type `T` and distance of type `D`.
"""
abstract type Variogram{T,D} end

"""
    isstationary(γ)

Check if variogram `γ` possesses the 2nd-order stationary property.
"""
isstationary(γ::Variogram) = isstationary(typeof(γ))

"""
    range(γ)

Return the range of the variogram `γ` when defined.
"""
Base.range(γ::Variogram) = γ.range

"""
    sill(γ)

Return the sill of the variogram `γ` when defined.
"""
sill(γ::Variogram) = γ.sill

"""
    nugget(γ)

Return the nugget of the variogram `γ` when defined.
"""
nugget(γ::Variogram) = γ.nugget

"""
    distance(γ)

Return the distance of the variogram `γ`.
"""
distance(γ::Variogram) = γ.distance

"""
    γ(u, v)

Evaluate the variogram at points `u` and `v`.
"""
(γ::Variogram)(u::Point, v::Point) =
  γ(evaluate(distance(γ), coordinates(u), coordinates(v)))

"""
    γ(U, v)

Evaluate the variogram at geometry `U` and point `v`.
"""
function (γ::Variogram)(U::Geometry, v::Point)
  us = _sample(γ, U)
  mean(γ(u, v) for u in us)
end

"""
    γ(u, V)

Evaluate the variogram at point `u` and geometry `V`.
"""
(γ::Variogram)(u::Point, V::Geometry) = γ(V, u)

"""
    γ(U, V)

Evaluate the variogram at geometries `U` and `V`.
"""
function (γ::Variogram)(U::Geometry, V::Geometry)
  us = _sample(γ, U)
  vs = _sample(γ, V)
  mean(γ(u, v) for u in us, v in vs)
end

"""
    result_type(γ, u, v)

Return result type of γ(u, v).
"""
result_type(γ::Variogram, u, v) = typeof(γ(u, v))

#------------------
# IMPLEMENTATIONS
#------------------

"""
    GaussianVariogram(sill=s, range=r, nugget=n, distance=d)

A Gaussian variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct GaussianVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.0
  sill::T   = 1.0
  nugget::T = 0.0
  distance::D = Euclidean()
end
function (γ::GaussianVariogram)(h)
  # add small eps to nugget
  # for numerical stability
  s = γ.sill
  r = γ.range
  n = γ.nugget + typeof(s)(1e-6)
  (s - n) * (1 - exp(-3(h/r)^2)) + (h > 0) * n
end
isstationary(::Type{<:GaussianVariogram}) = true

"""
    ExponentialVariogram(sill=s, range=r, nugget=n, distance=d)

An exponential variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct ExponentialVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.0
  sill::T   = 1.0
  nugget::T = 0.0
  distance::D = Euclidean()
end
function (γ::ExponentialVariogram)(h)
  s = γ.sill
  r = γ.range
  n = γ.nugget
  (s - n) * (1 - exp(-3(h/r))) + (h > 0) * n
end
isstationary(::Type{<:ExponentialVariogram}) = true

"""
    MaternVariogram(sill=s, range=r, nugget=n, order=ν, distance=d)

A Matérn variogram with sill `s`, range `r` and nugget `n`. The parameter
ν is the order of the Bessel function. Optionally, use a custom distance `d`.
"""
@with_kw struct MaternVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.0
  sill::T   = 1.0
  nugget::T = 0.0
  order::T  = 1.0
  distance::D = Euclidean()
end
function (γ::MaternVariogram)(h::T) where {T}
  s = γ.sill
  r = γ.range
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

"""
    SphericalVariogram(sill=s, range=r, nugget=n, distance=d)

A spherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SphericalVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.0
  sill::T   = 1.0
  nugget::T = 0.0
  distance::D = Euclidean()
end
function (γ::SphericalVariogram)(h::T) where {T}
  s = γ.sill
  r = γ.range
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

"""
    CubicVariogram(sill=s, range=r, nugget=n, distance=d)

A cubic variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct CubicVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.0
  sill::T   = 1.0
  nugget::T = 0.0
  distance::D = Euclidean()
end
function (γ::CubicVariogram)(h::T) where {T}
  s = γ.sill
  r = γ.range
  n = γ.nugget

  # constants
  c1 = T(35) / T(4)
  c2 = T(7)  / T(2)
  c3 = T(3)  / T(4)
  s1 = 7*(h/r)^2 - c1*(h/r)^3 + c2*(h/r)^5 - c3*(h/r)^7
  s2 = T(1)

  (h < r) * (s - n) * s1 +
  (h ≥ r) * (s - n) * s2 +
  (h > 0) * n
end
isstationary(::Type{<:CubicVariogram}) = true

"""
    PentasphericalVariogram(sill=s, range=r, nugget=n, distance=d)

A pentaspherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PentasphericalVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.0
  sill::T   = 1.0
  nugget::T = 0.0
  distance::D = Euclidean()
end
function (γ::PentasphericalVariogram)(h::T) where {T}
  s = γ.sill
  r = γ.range
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

"""
    PowerVariogram(scaling=s, exponent=a, nugget=n, distance=d)

A power variogram with scaling `s`, exponent `a` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PowerVariogram{T,D} <: Variogram{T,D}
  scaling::T  = 1.0
  nugget::T   = 0.0
  exponent::T = 1.0
  distance::D = Euclidean()
end
function (γ::PowerVariogram)(h)
  s = γ.scaling
  a = γ.exponent
  n = γ.nugget
  s * h^a + (h > 0) * n
end
isstationary(::Type{<:PowerVariogram}) = false

"""
    SineHoleVariogram(sill=s, range=r, nugget=n, distance=d)

A sine hole variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SineHoleVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.0
  sill::T   = 1.0
  nugget::T = 0.0
  distance::D = Euclidean()
end
function (γ::SineHoleVariogram)(h::T) where {T}
  s = γ.sill
  r = γ.range
  n = γ.nugget

  # shift lag by machine precision to
  # avoid explosion at the origin
  h′ = h + eps(T)
  c  = T(π)

  (s - n) * (1 - sin(c*h′/r)/(c*h′/r)) + (h′ > 0) * n
end
isstationary(::Type{<:SineHoleVariogram}) = true

"""
    NuggetEffect(n)
    NuggetEffect(nugget=n, distance=d)

A pure nugget effect variogram with nugget `n`.
Optionally use a custom distance `d`.
"""
@with_kw struct NuggetEffect{T,D} <: Variogram{T,D}
  nugget::T = 0.0
  distance::D = Euclidean()
end
NuggetEffect(n) = NuggetEffect(nugget=n)
function (γ::NuggetEffect)(h)
  (h > 0) * γ.nugget
end
Base.range(::NuggetEffect{T}) where {T} = zero(T)
sill(γ::NuggetEffect) = γ.nugget
isstationary(::Type{<:NuggetEffect}) = true
