# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    pairwise(γ, X)

Evaluate variogram `γ` between all n² pairs of columns in a
m-by-n matrix `X` efficiently.
"""
function pairwise(γ::Variogram, X::AbstractMatrix)
  m, n = size(X)
  Γ = Array{result_type(γ, X, X)}(n, n)
  for j=1:n
    xj = view(X, :, j)
    for i=j+1:n
      xi = view(X, :, i)
      @inbounds Γ[i,j] = γ(xi, xj)
    end
    @inbounds Γ[j,j] = γ(xj, xj)
    for i=1:j-1
      @inbounds Γ[i,j] = Γ[j,i] # leverage the symmetry
    end
  end

  Γ
end

"""
    pairwise(γ, domain, locations)

Evaluate variogram `γ` between all `locations` in `domain`.
"""
function pairwise(γ::Variogram, domain::AbstractDomain, locations::AbstractVector{Int})
  n = length(locations)
  x = coordinates(domain, locations[1])
  Γ = Array{result_type(γ, x, x)}(n, n)
  for j=1:n
    xj = coordinates(domain, locations[j])
    for i=j+1:n
      xi = coordinates(domain, locations[i])
      @inbounds Γ[i,j] = γ(xi, xj)
    end
    @inbounds Γ[j,j] = γ(xj, xj)
    for i=1:j-1
      @inbounds Γ[i,j] = Γ[j,i] # leverage the symmetry
    end
  end

  Γ
end

"""
    pairwise(γ, domain, locations₁, locations₂)

Evaluate variogram `γ` between `locations₁` and `locations₂` in `domain`.
"""
function pairwise(γ::Variogram, domain::AbstractDomain,
                  locations₁::AbstractVector{Int},
                  locations₂::AbstractVector{Int})
  m = length(locations₁)
  n = length(locations₂)
  x₁ = coordinates(domain, locations₁[1])
  x₂ = coordinates(domain, locations₂[1])
  Γ = Array{result_type(γ, x₁, x₂)}(m, n)
  for j=1:n
    xj = coordinates(domain, locations₂[j])
    for i=1:m
      xi = coordinates(domain, locations₁[i])
      @inbounds Γ[i,j] = γ(xi, xj)
    end
  end

  Γ
end
