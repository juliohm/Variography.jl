# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Covariance

Parent type of all covariance functions (e.g. Gaussian covariance).
"""
abstract type Covariance end

"""
    cov(x₁, x₂)

Evaluate the covariance at objects `x₁` and `x₁`.
"""
(cov::Covariance)(x₁, x₂) = sill(cov.γ) - cov.γ(x₁, x₂)

"""
    pairwise(cov, domain)
    
Evaluate covariance `cov` between all elements in the `domain`.
    
    pairwise(cov, domain₁, domain₂)

Evaluate covariance `cov` between all elements of `domain₁` and `domain₂`.
"""
pairwise(cov::Covariance, args...) = sill(cov.γ) .- pairwise(cov.γ, args...)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, cov::T) where {T<:Covariance}
  name = string(nameof(T))
  _showcompact(io, name, cov.γ)
end

function Base.show(io::IO, ::MIME"text/plain", cov::T) where {T<:Covariance}
  name = string(nameof(T))
  _showfull(io, name, cov.γ)
end

# heper macro to define covariances
macro defcov(CovType, VarioType)
  docstring = """
      $(CovType)(args..., kwargs...)

  A covariance function derived from the corresponding variogram function.

  Please see [`$(VarioType)`](@ref) for available parameters.
  """
  expr = quote
    @doc $docstring struct $CovType{V<:$VarioType} <: Covariance
      γ::V
    end

    $CovType(args...; kwargs...) = $CovType($VarioType(args...; kwargs...))
  end
  esc(expr)
end

@defcov CircularCovariance CircularVariogram

@defcov CubicCovariance CubicVariogram

@defcov ExponentialCovariance ExponentialVariogram

@defcov GaussianCovariance GaussianVariogram

@defcov MaternCovariance MaternVariogram

@defcov PentasphericalCovariance PentasphericalVariogram

@defcov SineHoleCovariance SineHoleVariogram

@defcov SphericalCovariance SphericalVariogram
