# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# recommended maximum number of points inside a geometry
# Journel, A. & Huijbregts, Ch.J. Mining Geostatistics page 97
function _maxpoints(g::Geometry)
  paramdim(g) == 1 && return 10
  paramdim(g) == 2 && return 6*6
  paramdim(g) == 3 && return 4*4*4
end

function (γ::Variogram)(U::Geometry, v::Point)
  μ = measure(U)
  d = paramdim(U)
  n = _maxpoints(U)
  α = (μ / n) ^ (1 / d)
  us = sample(U, MinDistanceSampling(α))
  mean(γ(u, v) for u in us)
end

(γ::Variogram)(u::Point, V::Geometry) = γ(V, u)

function (γ::Variogram)(U::Geometry, V::Geometry)
  μᵤ = measure(U)
  dᵤ = paramdim(U)
  nᵤ = _maxpoints(U)
  μᵥ = measure(V)
  dᵥ = paramdim(V)
  nᵥ = _maxpoints(V)
  αᵤ = (μᵤ / nᵤ) ^ (1 / dᵤ)
  αᵥ = (μᵥ / nᵥ) ^ (1 / dᵥ)
  us = sample(U, MinDistanceSampling(αᵤ))
  vs = sample(V, MinDistanceSampling(αᵥ))
  mean(γ(u, v) for u in us, v in vs)
end
