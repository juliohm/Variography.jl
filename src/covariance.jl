# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Covariance

Parent type of all covariance models (e.g. Gaussian covariance).
"""
abstract type Covariance end

"""
    cov(x₁, x₂)

Evaluate the covariance at objects `x₁` and `x₁`.
"""
(cov::Covariance)(x₁, x₂) = sill(cov.γ) - cov.γ(x₁, x₂)

# heper macro to define covariances
macro defcov(CovType, VarioType)
  expr = quote
    Base.@__doc__ struct $CovType{V<:$VarioType} <: Covariance
      γ::V
    end

    $CovType(args...; kwargs...) = $CovType($VarioType(args...; kwargs...))
  end
  esc(expr)
end

"""
    CircularCovariance(; range=r, sill=s, nugget=n)
    CircularCovariance(ball; sill=s, nugget=n)

A circular covariance with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
@defcov CircularCovariance CircularVariogram 

"""
    CubicCovariance(range=r, sill=s, nugget=n)
    CubicCovariance(ball; sill=s, nugget=n)

A cubic covariance with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
@defcov CubicCovariance CubicVariogram 

"""
    ExponentialCovariance(range=r, sill=s, nugget=n)
    ExponentialCovariance(ball; sill=s, nugget=n)

An exponential covariance with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
@defcov ExponentialCovariance ExponentialVariogram 

"""
    GaussianCovariance(range=r, sill=s, nugget=n)
    GaussianCovariance(ball; sill=s, nugget=n)

A Gaussian covariance with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
@defcov GaussianCovariance GaussianVariogram 

"""
    MaternCovariance(range=r, sill=s, nugget=n, order=ν)
    MaternCovariance(ball; sill=s, nugget=n, order=ν)

A Matérn covariance with range `r`, sill `s` and nugget `n`.
The parameter `ν` is the order of the Bessel function.
Optionally, use a custom metric `ball`.
"""
@defcov MaternCovariance MaternVariogram 

"""
    NuggetEffect(nugget=n)

A pure nugget effect covariance with nugget `n`.
"""
@defcov NuggetCovariance NuggetEffect 

"""
    PentasphericalCovariance(range=r, sill=s, nugget=n)
    PentasphericalCovariance(ball; sill=s, nugget=n)

A pentaspherical covariance with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
@defcov PentasphericalCovariance PentasphericalVariogram 

"""
    SineHoleCovariance(range=r, sill=s, nugget=n)
    SineHoleCovariance(ball; sill=s, nugget=n)

A sine hole covariance with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
@defcov SineHoleCovariance SineHoleVariogram 

"""
    SphericalCovariance(range=r, sill=s, nugget=n)
    SphericalCovariance(ball; sill=s, nugget=n)

A spherical covariance with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
@defcov SphericalCovariance SphericalVariogram 
