# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function (γ::Variogram)(U::Geometry, v::Point)
  us = _reg_sample(U)
  mean(γ(u, v) for u in us)
end

(γ::Variogram)(u::Point, V::Geometry) = γ(V, u)

function (γ::Variogram)(U::Geometry, V::Geometry)
  us = _reg_sample(U)
  vs = _reg_sample(V)
  mean(γ(u, v) for u in us, v in vs)
end

# helper function to sample points within a geometry
function _reg_sample(g::Geometry)
  α = _reg_spacing(g)
  sample(g, MinDistanceSampling(α))
end

# helper function to compute the spacing of points within
# a given geometry based on a maximum number of points
function _reg_spacing(g::Geometry)
  μ = measure(g)
  d = paramdim(g)
  n = _reg_maxpoints(g)
  (μ / n) ^ (1 / d)
end

# recommended maximum number of points inside a geometry
# Journel, A. & Huijbregts, Ch.J. Mining Geostatistics page 97
function _reg_maxpoints(g::Geometry)
  paramdim(g) == 1 && return 10
  paramdim(g) == 2 && return 6*6
  paramdim(g) == 3 && return 4*4*4
end

# --------------
# SPECIAL CASES
# --------------

function _reg_sample(s::Segment)
  s′ = Segment(s(.05), s(.95))
  sample(s′, RegularSampling(9))
end

function _reg_sample(q::Quadrangle)
  q′ = Quadrangle(q(.05,.05), q(.95,.05), q(.95,.95), q(.05,.95))
  sample(q′, RegularSampling(3, 3))
end

function _reg_sample(h::Hexahedron)
  h′ = Hexahedron(h(.05,.05,.05), h(.95,.05,.05), h(.95,.95,.05), h(.05,.95,.05),
                  h(.05,.05,.95), h(.95,.05,.95), h(.95,.95,.95), h(.05,.95,.95))
  sample(h′, RegularSampling(3, 3, 3))
end
