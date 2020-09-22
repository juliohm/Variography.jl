# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Variogram{T,D}

A theoretical variogram model (e.g. Gaussian variogram)
with parameters of type `T` and distance of type `D`.
"""
abstract type Variogram{T,D} end

"""
    result_type(γ, x₁, x₂)

Return result type of γ(x₁, x₂).
"""
result_type(γ::Variogram, x₁::AbstractArray, x₂::AbstractArray) = typeof(γ(x₁, x₂))

"""
    isstationary(γ)
    isstationary(V)

Check if variogram `γ` or variogram type `V` possesses
the 2nd-order stationary property.
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
    γ(x, y)

Evaluate the variogram at points `x` and `y`.
"""
(γ::Variogram)(x, y) = γ(evaluate(γ.distance, x, y))

#------------------
# IMPLEMENTATIONS
#------------------
"""
    GaussianVariogram(sill=s, range=r, nugget=n, distance=d)

A Gaussian variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct GaussianVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 1e-6 # positive nugget for numerical stability
  distance::D = Euclidean()
end
(γ::GaussianVariogram)(h) = (γ.sill - γ.nugget) * (1 - exp(-3(h/γ.range)^2)) + (h > 0) * γ.nugget
isstationary(::Type{<:GaussianVariogram}) = true

"""
    ExponentialVariogram(sill=s, range=r, nugget=n, distance=d)

An exponential variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct ExponentialVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::ExponentialVariogram)(h) = (γ.sill - γ.nugget) * (1 - exp(-3(h/γ.range))) + (h > 0) * γ.nugget
isstationary(::Type{<:ExponentialVariogram}) = true

"""
    MaternVariogram(sill=s, range=r, nugget=n, order=ν, distance=d)

A Matérn variogram with sill `s`, range `r` and nugget `n`. The parameter
ν is the order of the Bessel function. Optionally, use a custom distance `d`.
"""
@with_kw struct MaternVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  order::T  = 1.
  distance::D = Euclidean()
end
(γ::MaternVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget
  ν = γ.order

  # shift lag by machine precision to
  # avoid explosion at the origin
  h2 = h + eps(eltype(h))
  h3 = sqrt(2.0ν)h2/r

  (s - n) * (1 - 2.0^(1 - ν)/gamma(ν) * h3^ν * besselk(ν, h3)) + (h > 0) * n
end
isstationary(::Type{<:MaternVariogram}) = true

"""
    SphericalVariogram(sill=s, range=r, nugget=n, distance=d)

A spherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SphericalVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::SphericalVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget

  (h < r) * (s - n) * (1.5(h/r) - 0.5(h/r)^3) + (h ≥ r) * (s - n) + (h > 0) * n
end
isstationary(::Type{<:SphericalVariogram}) = true

"""
    CubicVariogram(sill=s, range=r, nugget=n, distance=d)

A cubic variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct CubicVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::CubicVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget

  (h < r) * (s - n) * (7*(h/r)^2 - (35/4)*(h/r)^3 + (7/2)*(h/r)^5 - (3/4)*(h/r)^7) +
  (h ≥ r) * (s - n) + (h > 0) * n
end
isstationary(::Type{<:CubicVariogram}) = true

"""
    PentasphericalVariogram(sill=s, range=r, nugget=n, distance=d)

A pentaspherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PentasphericalVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::PentasphericalVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget

  (h < r) * (s - n) * ((15/8)*(h/r) - (5/4)*(h/r)^3 + (3/8)*(h/r)^5) +
  (h ≥ r) * (s - n) + (h > 0) * n
end
isstationary(::Type{<:PentasphericalVariogram}) = true

"""
    PowerVariogram(scaling=s, exponent=a, nugget=n, distance=d)

A power variogram with scaling `s`, exponent `a` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PowerVariogram{T,D} <: Variogram{T,D}
  scaling::T  = 1.
  nugget::T   = 0.
  exponent::T = 1.
  distance::D = Euclidean()
end
(γ::PowerVariogram)(h) = γ.scaling*h^γ.exponent + (h > 0) * γ.nugget
isstationary(::Type{<:PowerVariogram}) = false

"""
    SineHoleVariogram(sill=s, range=r, nugget=n, distance=d)

A sine hole variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SineHoleVariogram{T,D} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::SineHoleVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget

  # shift lag by machine precision to
  # avoid explosion at the origin
  h = h + eps(eltype(h))

  (s - n) * (1 - sin(π*h/r)/(π*h/r)) + (h > 0) * n
end
isstationary(::Type{<:SineHoleVariogram}) = true

#----------------------
# VARIOGRAM OPERATIONS
#----------------------
"""
    SumVariogram(γ₁, γ₂)

Internal representation of the sum `γ₁ + γ₂`.
"""
struct SumVariogram{G₁,G₂} <: Variogram{Number,Metric}
  γ₁::G₁
  γ₂::G₂
end
(γ::SumVariogram)(h) = γ.γ₁(h) + γ.γ₂(h)
(γ::SumVariogram)(x, y) = γ.γ₁(x, y) + γ.γ₂(x, y)
Base.range(γ::SumVariogram) = max(range(γ.γ₁), range(γ.γ₂))
sill(γ::SumVariogram) = sill(γ.γ₁) + sill(γ.γ₂)
nugget(γ::SumVariogram) = nugget(γ.γ₁) + nugget(γ.γ₂)
isstationary(γ::SumVariogram) = isstationary(γ.γ₁) && isstationary(γ.γ₂)

"""
    ScaledVariogram(c, γ)

Internal representation of the scalar multiplication `cγ`.
"""
struct ScaledVariogram{T,D,G<:Variogram{T,D},V} <: Variogram{T,D}
  c::V
  γ::G
end
(γ::ScaledVariogram)(h) = γ.c * γ.γ(h)
(γ::ScaledVariogram)(x, y) = γ.c * γ.γ(x, y)
Base.range(γ::ScaledVariogram) = range(γ.γ)
sill(γ::ScaledVariogram) = γ.c * sill(γ.γ)
nugget(γ::ScaledVariogram) = γ.c * nugget(γ.γ)
isstationary(γ::ScaledVariogram) = isstationary(γ.γ)

"""
    γ₁ + γ₂

Return sum of variograms `γ₁` and `γ₂`.
"""
+(γ₁::Variogram, γ₂::Variogram) = SumVariogram(γ₁, γ₂)

"""
    cγ

Return the multiplication of the scalar `c` with the variogram `γ`.
"""
*(c::Number, γ::Variogram) = ScaledVariogram(c, γ)

# collect all stationary models for other parts of the codebase
const OPERATIONS = [SumVariogram, ScaledVariogram]
const STATIONARY = filter(isstationary, setdiff(subtypes(Variogram), OPERATIONS))
