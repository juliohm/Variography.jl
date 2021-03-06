# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    pairwise(γ, X)

Evaluate variogram `γ` between all n² pairs of columns in a
m-by-n matrix `X` efficiently.
"""
function pairwise(γ::Variogram, X::AbstractMatrix)
  m, n = size(X)
  R = result_type(γ, X, X)
  Γ = Array{R}(undef, n, n)
  pairwise!(Γ, γ, X)

  Γ
end

"""
    pairwise!(Γ, γ, X)

Non-allocating pairwise evaluation.
"""
function pairwise!(Γ, γ::Variogram, X::AbstractMatrix)
  m, n = size(X)
  @inbounds for j=1:n
    xj = view(X, :, j)
    for i=j+1:n
      xi = view(X, :, i)
      Γ[i,j] = γ(xi, xj)
    end
    Γ[j,j] = γ(xj, xj)
    for i=1:j-1
      Γ[i,j] = Γ[j,i] # leverage the symmetry
    end
  end
end

"""
    pairwise(γ, domain, inds)

Evaluate variogram `γ` between all `inds` in `domain`.
"""
function pairwise(γ::Variogram, domain,
                  inds::AbstractVector{Int})
  Dim = embeddim(domain)
  T = coordtype(domain)
  x = rand(SVector{Dim,T})
  n = length(inds)
  R = result_type(γ, x, x)
  Γ = Matrix{R}(undef, n, n)
  @inbounds for j=1:n
    xj = coordinates(centroid(domain, inds[j]))
    for i=j+1:n
      xi = coordinates(centroid(domain, inds[i]))
      Γ[i,j] = γ(xi, xj)
    end
    Γ[j,j] = γ(xj, xj)
    for i=1:j-1
      Γ[i,j] = Γ[j,i] # leverage the symmetry
    end
  end

  Γ
end

"""
    pairwise(γ, domain, inds₁, inds₂)

Evaluate variogram `γ` between `inds₁` and `inds₂` in `domain`.
"""
function pairwise(γ::Variogram, domain,
                  inds₁::AbstractVector{Int},
                  inds₂::AbstractVector{Int})
  Dim = embeddim(domain)
  T = coordtype(domain)
  x = rand(SVector{Dim,T})
  m = length(inds₁)
  n = length(inds₂)
  R = result_type(γ, x, x)
  Γ = Array{R}(undef, m, n)
  @inbounds for j=1:n
    xj = coordinates(centroid(domain, inds₂[j]))
    for i=1:m
      xi = coordinates(centroid(domain, inds₁[i]))
      Γ[i,j] = γ(xi, xj)
    end
  end

  Γ
end
