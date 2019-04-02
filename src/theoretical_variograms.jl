# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    Variogram

A variogram model (e.g. Gaussian variogram).
"""
abstract type Variogram{T,D} end

"""
    isstationary(γ)

Check if variogram `γ` possesses the 2nd-order stationary property.
"""
isstationary(::Variogram) = false

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
    param_type(γ)

Return the parameter (e.g. sill, range) type of the variogram.
"""
param_type(::Variogram{T,D}) where {T<:Real,D<:Metric} = T

"""
    result_type(γ, x₁, x₂)

Return result type of γ(x₁, x₂).
"""
result_type(γ::Variogram, x₁::AbstractArray, x₂::AbstractArray) =
  promote_type(param_type(γ), Distances.result_type(γ.distance, x₁, x₂))

#------------------
# IMPLEMENTATIONS
#------------------
"""
    GaussianVariogram(sill=s, range=r, nugget=n, distance=d)

A Gaussian variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct GaussianVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 1e-6 # positive nugget for numerical stability
  distance::D = Euclidean()
end
(γ::GaussianVariogram)(h) = (γ.sill - γ.nugget) * (1 - exp(-3(h/γ.range)^2)) + (h > 0) * γ.nugget
(γ::GaussianVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::GaussianVariogram) = true

"""
    ExponentialVariogram(sill=s, range=r, nugget=n, distance=d)

An exponential variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct ExponentialVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::ExponentialVariogram)(h) = (γ.sill - γ.nugget) * (1 - exp(-3(h/γ.range))) + (h > 0) * γ.nugget
(γ::ExponentialVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::ExponentialVariogram) = true

"""
    MaternVariogram(sill=s, range=r, nugget=n, order=ν, distance=d)

A Matérn variogram with sill `s`, range `r` and nugget `n`. The parameter
ν is the order of the Bessel function. Optionally, use a custom distance `d`.
"""
@with_kw struct MaternVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
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
(γ::MaternVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::MaternVariogram) = true

"""
    SphericalVariogram(sill=s, range=r, nugget=n, distance=d)

A spherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SphericalVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
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
(γ::SphericalVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::SphericalVariogram) = true

"""
    CubicVariogram(sill=s, range=r, nugget=n, distance=d)

A cubic variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct CubicVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
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
(γ::CubicVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::CubicVariogram) = true

"""
    PentasphericalVariogram(sill=s, range=r, nugget=n, distance=d)

A pentaspherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PentasphericalVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
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
(γ::PentasphericalVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::PentasphericalVariogram) = true

"""
    PowerVariogram(scaling=s, exponent=a, nugget=n, distance=d)

A power variogram with scaling `s`, exponent `a` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PowerVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
  scaling::T  = 1.
  nugget::T   = 0.
  exponent::T = 1.
  distance::D = Euclidean()
end
(γ::PowerVariogram)(h) = @. γ.scaling*h^γ.exponent + (h > 0) * γ.nugget
(γ::PowerVariogram)(x, y) = γ(evaluate(γ.distance, x, y))

"""
    SineHoleVariogram(sill=s, range=r, nugget=n, distance=d)

A sine hole variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SineHoleVariogram{T<:Real,D<:Metric} <: Variogram{T,D}
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
(γ::SineHoleVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::SineHoleVariogram) = true

#------------------------
# COMPOSITIVE VARIOGRAMS
#------------------------
"""
    CompositeVariogram(γ₁, γ₂, ..., γₙ)

A composite (additive) model of variograms γ(h) = γ₁(h) + γ₂(h) + ⋯ + γₙ(h).
"""
struct CompositeVariogram <: Variogram{Real,Metric}
  γs::Vector{Variogram}
  CompositeVariogram(g, gs...) = new([g, gs...])
end
(c::CompositeVariogram)(h) = sum(γ(h) for γ in c.γs)
(c::CompositeVariogram)(x, y) = sum(γ(x,y) for γ in c.γs)
isstationary(c::CompositeVariogram) = all(isstationary(γ) for γ in c.γs)
sill(c::CompositeVariogram) = sum(sill(γ) for γ in c.γs)
param_type(c::CompositeVariogram) = promote_type([param_type(γ) for γ in c.γs]...)

"""
    γ₁ + γ₂

Return compositive (additive) model of variograms `γ₁` and `γ₂`.
"""
+(γ₁::Variogram, γ₂::Variogram) = CompositeVariogram(γ₁, γ₂)
+(c::CompositeVariogram, γ::Variogram) = CompositeVariogram(c.γs..., γ)
+(γ::Variogram, c::CompositeVariogram) = CompositeVariogram(γ, c.γs...)
+(c₁::CompositeVariogram, c₂::CompositeVariogram) = CompositeVariogram(c₁.γs..., c₂.γs...)
