# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function (γ::Variogram)(U::Geometry, v::Point)
  us = _reg_sample(γ, U)
  mean(γ(u, v) for u in us)
end

(γ::Variogram)(u::Point, V::Geometry) = γ(V, u)

function (γ::Variogram)(U::Geometry, V::Geometry)
  us = _reg_sample(γ, U)
  vs = _reg_sample(γ, V)
  mean(γ(u, v) for u in us, v in vs)
end

# helper function to sample points within a geometry
# given a variogram model with well-defined range
function _reg_sample(γ::Variogram, V::Geometry)
  α = _reg_spacing(γ, V)
  sample(V, MinDistanceSampling(α))
end

function _reg_spacing(γ::Variogram, V::Geometry)
  s = sides(boundingbox(V))
  l = filter(>(0), s)
  min(range(γ), minimum(l)) / 3
end

function _reg_dims(γ::Variogram, V::Geometry)
  s = sides(boundingbox(V))
  α = _reg_spacing(γ, V)
  n = ceil.(Int, s ./ α)
  replace(n, 0 => 1)
end

# --------------
# SPECIAL CASES
# --------------

_reg_sample(::Variogram, p::Point) = [p]

function _reg_sample(γ::Variogram, s::Segment)
  s′ = Segment(s(.05), s(.95))
  n  = _reg_dims(γ, s)
  sample(s′, RegularSampling(n...))
end

function _reg_sample(γ::Variogram, q::Quadrangle)
  q′ = Quadrangle(q(.05,.05), q(.95,.05), q(.95,.95), q(.05,.95))
  n  = _reg_dims(γ, q)
  sample(q′, RegularSampling(n...))
end

function _reg_sample(γ::Variogram, h::Hexahedron)
  h′ = Hexahedron(h(.05,.05,.05), h(.95,.05,.05), h(.95,.95,.05), h(.05,.95,.05),
                  h(.05,.05,.95), h(.95,.05,.95), h(.95,.95,.95), h(.05,.95,.95))
  n  = _reg_dims(γ, h)
  sample(h′, RegularSampling(n...))
end