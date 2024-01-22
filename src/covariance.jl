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

# heper macro to define covariances
macro defcov(CovType, VarioType)
  docstring = """
      $(CovType)(args..., kwargs...)

  A covariance function derived from the corresponding variogram function.

  Please see [`$(VarioType)`](@ref) for available parameters.
  """
  expr = quote
    @doc $docstring
    struct $CovType{V<:$VarioType} <: Covariance
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

@defcov NuggetCovariance NuggetEffect 

@defcov PentasphericalCovariance PentasphericalVariogram 

@defcov SineHoleCovariance SineHoleVariogram 

@defcov SphericalCovariance SphericalVariogram 
