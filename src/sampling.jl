# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

_sample(::Variogram, p::Point) = [p]

function _sample(γ::Variogram, g::Geometry)
  α = _spacing(γ, g)
  sample(g, MinDistanceSampling(α))
end

function _sample(γ::Variogram, s::Segment)
  s′ = Segment(s(.05), s(.95))
  n  = _dims(γ, s)
  sample(s′, RegularSampling(n...))
end

function _sample(γ::Variogram, q::Quadrangle)
  q′ = Quadrangle(q(.05,.05), q(.95,.05), q(.95,.95), q(.05,.95))
  n  = _dims(γ, q)
  sample(q′, RegularSampling(n...))
end

function _sample(γ::Variogram, h::Hexahedron)
  h′ = Hexahedron(h(.05,.05,.05), h(.95,.05,.05), h(.95,.95,.05), h(.05,.95,.05),
                  h(.05,.05,.95), h(.95,.05,.95), h(.95,.95,.95), h(.05,.95,.95))
  n  = _dims(γ, h)
  sample(h′, RegularSampling(n...))
end

function _spacing(γ::Variogram, g::Geometry)
  s = sides(boundingbox(g))
  l = filter(>(0), s)
  min(range(γ), minimum(l)) / 3
end

function _dims(γ::Variogram, g::Geometry)
  s = sides(boundingbox(g))
  α = _spacing(γ, g)
  n = ceil.(Int, s ./ α)
  replace(n, 0 => 1)
end
