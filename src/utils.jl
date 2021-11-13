# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    spheredir(θ, φ)

Return the 3D direction given polar angle `θ` and
azimuthal angle `φ` in degrees according to the ISO
convention.
"""
function spheredir(theta, phi)
  θ, φ = deg2rad(theta), deg2rad(phi)
  Vec(sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ))
end

"""
    planebasis(normal)

Return 2D basis vectors in the plane with given 3D `normal`.
"""
function planebasis(normal::Vec{3,T}) where {T}
  # normalize input
  n = normal ./ sqrt(sum(normal[i]^2 for i in 1:3))

  # find last non-zero component
  idx = -1
  for (i, c) in enumerate(reverse(n))
    if c != 0
      idx = length(n) - i + 1
      break
    end
  end

  @assert idx > 0 "invalid normal vector"

  # first basis vector (perturb and subtract projection)
  u = ntuple(i -> i == idx%3 + 1 ? n[i] + one(T) : n[i], 3)
  l = sum(u[i]*n[i] for i in 1:3)
  u = ntuple(i -> u[i] - l*n[i], 3)

  # second basis vector (cross product)
  nx, ny, nz = n
  ux, uy, uz = u
  v = (ny*uz - nz*uz, nz*ux - nx*uz, nx*uy - ny*ux)

  # normalize output
  u = u ./ sqrt(sum(u[i]^2 for i in 1:3))
  v = v ./ sqrt(sum(v[i]^2 for i in 1:3))

  Vec(u), Vec(v)
end

planebasis(normal::NTuple{3,T}) where {T} = planebasis(Vec(normal))
