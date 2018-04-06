# ------------------------------------------------------------------
# Copyright (c) 2015, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    Variogram

A variogram model (e.g. Gaussian variogram).
"""
abstract type Variogram end

"""
    isstationary(γ)

Check if variogram `γ` possesses the 2nd-order stationary property.
"""
isstationary(::Variogram) = false

"""
    pairwise(γ, X)

Evaluate variogram `γ` between all n² pairs of columns in a
m-by-n matrix `X` efficiently.
"""
function pairwise(γ::Variogram, X::AbstractMatrix)
  m, n = size(X)
  Γ = Array{Float64}(n, n)
  for j=1:n
    xj = view(X, :, j)
    for i=j+1:n
      xi = view(X, :, i)
      @inbounds Γ[i,j] = γ(xi, xj)
    end
    @inbounds Γ[j,j] = γ(xj, xj)
    for i=1:j-1
      @inbounds Γ[i,j] = Γ[j,i] # leverage the symmetry
    end
  end

  Γ
end

#------------------
# IMPLEMENTATIONS
#------------------
"""
    GaussianVariogram(sill=s, range=r, nugget=n, distance=d)

A Gaussian variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct GaussianVariogram{T<:Real,D<:Metric} <: Variogram
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::GaussianVariogram)(h) = (γ.sill - γ.nugget) * (1 - exp.(-3(h/γ.range).^2)) + γ.nugget
(γ::GaussianVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::GaussianVariogram) = true

"""
    ExponentialVariogram(sill=s, range=r, nugget=n, distance=d)

An exponential variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct ExponentialVariogram{T<:Real,D<:Metric} <: Variogram
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::ExponentialVariogram)(h) = (γ.sill - γ.nugget) * (1 - exp.(-3(h/γ.range))) + γ.nugget
(γ::ExponentialVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::ExponentialVariogram) = true

"""
    MaternVariogram(sill=s, range=r, nugget=n, order=ν, distance=d)

A Matérn variogram with sill `s`, range `r` and nugget `n`. The parameter
ν is the order of the Bessel function. Optionally, use a custom distance `d`.
"""
@with_kw struct MaternVariogram{T<:Real,D<:Metric} <: Variogram
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
  h3 = sqrt.(2.0ν)h2/r

  (s - n) * (1 - 2.0^(1 - ν)/gamma(ν) * h3.^ν .* besselk.(ν, h3)) + n
end
(γ::MaternVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::MaternVariogram) = true

"""
    SphericalVariogram(sill=s, range=r, nugget=n, distance=d)

A spherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SphericalVariogram{T<:Real,D<:Metric} <: Variogram
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::SphericalVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget

  (h .< r) .* (s - n) .* (1.5(h/r) - 0.5(h/r).^3) + (h .≥ r) .* (s - n) + n
end
(γ::SphericalVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::SphericalVariogram) = true

"""
    CubicVariogram(sill=s, range=r, nugget=n, distance=d)

A cubic variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct CubicVariogram{T<:Real,D<:Metric} <: Variogram
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::CubicVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget

  (h .< r) .* (s - n) .* (7*(h/r).^2 - (35/4)*(h/r).^3 + (7/2)*(h/r).^5 - (3/4)*(h/r).^7) +
  (h .≥ r) .* (s - n) + n
end
(γ::CubicVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::CubicVariogram) = true

"""
    PentasphericalVariogram

A pentaspherical variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PentasphericalVariogram{T<:Real,D<:Metric} <: Variogram
  range::T  = 1.
  sill::T   = 1.
  nugget::T = 0.
  distance::D = Euclidean()
end
(γ::PentasphericalVariogram)(h) = begin
  s = γ.sill
  r = γ.range
  n = γ.nugget

  (h .< r) .* (s - n) .* ((15/8)*(h/r) - (5/4)*(h/r).^3 + (3/8)*(h/r).^5) +
  (h .≥ r) .* (s - n) + n
end
(γ::PentasphericalVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::PentasphericalVariogram) = true

"""
    PowerVariogram(scaling=s, exponent=a, nugget=n, distance=d)

A power variogram with scaling `s`, exponent `a` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct PowerVariogram{T<:Real,D<:Metric} <: Variogram
  scaling::T  = 1.
  nugget::T   = 0.
  exponent::T = 1.
  distance::D = Euclidean()
end
(γ::PowerVariogram)(h) = γ.scaling*h.^γ.exponent + γ.nugget
(γ::PowerVariogram)(x, y) = γ(evaluate(γ.distance, x, y))

"""
    SineHoleVariogram(sill=s, range=r, nugget=n, distance=d)

A sine hole variogram with sill `s`, range `r` and nugget `n`.
Optionally, use a custom distance `d`.
"""
@with_kw struct SineHoleVariogram{T<:Real,D<:Metric} <: Variogram
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

  (s - n) * (1 - sin.(π*h/r)./(π*h/r)) + n
end
(γ::SineHoleVariogram)(x, y) = γ(evaluate(γ.distance, x, y))
isstationary(::SineHoleVariogram) = true

"""
    CompositeVariogram(γ₁, γ₂, ..., γₙ)

A composite (additive) model of variograms γ(h) = γ₁(h) + γ₂(h) + ⋯ + γₙ(h).
"""
struct CompositeVariogram <: Variogram
  γs::Vector{Variogram}
  CompositeVariogram(g, gs...) = new([g, gs...])
end
(c::CompositeVariogram)(h) = sum(γ(h) for γ in c.γs)
(c::CompositeVariogram)(x, y) = sum(γ(x,y) for γ in c.γs)
