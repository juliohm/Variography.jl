# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

variosample(::Variogram, p::Point) = [p]

function variosample(γ::Variogram, g::Geometry)
  α = _spacing(γ, g)
  rng = MersenneTwister(123)
  sample(rng, g, MinDistanceSampling(α))
end

function variosample(γ::Variogram, s::Segment)
  s′ = Segment(s(0.05), s(0.95))
  n = _dims(γ, s)
  sample(s′, RegularSampling(n...))
end

function variosample(γ::Variogram, q::Quadrangle)
  q′ = Quadrangle(q(0.05, 0.05), q(0.95, 0.05), q(0.95, 0.95), q(0.05, 0.95))
  n = _dims(γ, q)
  sample(q′, RegularSampling(n...))
end

function variosample(γ::Variogram, h::Hexahedron)
  h′ = Hexahedron(
    h(0.05, 0.05, 0.05),
    h(0.95, 0.05, 0.05),
    h(0.95, 0.95, 0.05),
    h(0.05, 0.95, 0.05),
    h(0.05, 0.05, 0.95),
    h(0.95, 0.05, 0.95),
    h(0.95, 0.95, 0.95),
    h(0.05, 0.95, 0.95)
  )
  n = _dims(γ, h)
  sample(h′, RegularSampling(n...))
end

function _spacing(γ::Variogram, g::Geometry)
  r = range(γ)
  s = sides(boundingbox(g))
  l = minimum(filter(>(0), s))
  r > 0 ? min(r, l) / 3 : l / 3
end

function _dims(γ::Variogram, g::Geometry)
  s = sides(boundingbox(g))
  α = _spacing(γ, g)
  n = ceil.(Int, s ./ α)
  ntuple(i -> iszero(n[i]) ? 1 : n[i], length(n))
end
