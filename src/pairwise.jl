# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
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
    pairwise(γ, domain, locations)

Evaluate variogram `γ` between all `locations` in `domain`.
"""
function pairwise(γ::Variogram, domain::AbstractDomain{T,N},
                  locations::AbstractVector{Int}) where {T<:Real,N}
  xi = MVector{N,T}(undef)
  xj = MVector{N,T}(undef)
  n = length(locations)
  R = result_type(γ, xi, xj)
  Γ = Matrix{R}(undef, n, n)
  @inbounds for j=1:n
    coordinates!(xj, domain, locations[j])
    for i=j+1:n
      coordinates!(xi, domain, locations[i])
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
    pairwise(γ, domain, locations₁, locations₂)

Evaluate variogram `γ` between `locations₁` and `locations₂` in `domain`.
"""
function pairwise(γ::Variogram, domain::AbstractDomain{T,N},
                  locations₁::AbstractVector{Int},
                  locations₂::AbstractVector{Int}) where {T<:Real,N}
  xi = MVector{N,T}(undef)
  xj = MVector{N,T}(undef)
  m = length(locations₁)
  n = length(locations₂)
  R = result_type(γ, xi, xj)
  Γ = Array{R}(undef, m, n)
  @inbounds for j=1:n
    coordinates!(xj, domain, locations₂[j])
    for i=1:m
      coordinates!(xi, domain, locations₁[i])
      Γ[i,j] = γ(xi, xj)
    end
  end

  Γ
end
