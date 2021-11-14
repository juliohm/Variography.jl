# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Variogram

A theoretical variogram model (e.g. Gaussian variogram).
"""
abstract type Variogram end

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
    range(γ)

Return the maximum range of the variogram `γ` when defined.
"""
Base.range(γ::Variogram) = maximum(radii(γ.ball))

"""
    γ(u, v)

Evaluate the variogram at points `u` and `v`.
"""
function (γ::Variogram)(u::Point, v::Point)
  d = metric(γ.ball)
  x = coordinates(u)
  y = coordinates(v)
  h = evaluate(d, x, y)
  γ(h)
end

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

"""
    isstationary(γ)

Check if variogram `γ` possesses the 2nd-order stationary property.
"""
isstationary(γ::Variogram) = isstationary(typeof(γ))

#------------------
# IMPLEMENTATIONS
#------------------

include("theoretical/gaussian.jl")
include("theoretical/exponential.jl")
include("theoretical/spherical.jl")
include("theoretical/matern.jl")
include("theoretical/cubic.jl")
include("theoretical/pentaspherical.jl")
include("theoretical/sinehole.jl")
include("theoretical/power.jl")
include("theoretical/nugget.jl")