# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    pairwise(γ, domain)

Evaluate variogram `γ` between all elements in the `domain`.
"""
function pairwise(γ::Variogram, domain)
  u = first(domain)
  n = nelements(domain)
  R = result_type(γ, u, u)
  Γ = Matrix{R}(undef, n, n)
  pairwise!(Γ, γ, domain)
end

function pairwise!(Γ, γ::Variogram, domain)
  n = nelements(domain)
  @inbounds for j in 1:n
    vj = domain[j]
    for i in j+1:n
      ui = domain[i]
      Γ[i,j] = γ(ui, vj)
    end
    Γ[j,j] = γ(vj, vj)
    for i in 1:j-1
      Γ[i,j] = Γ[j,i] # leverage the symmetry
    end
  end
  Γ
end

"""
    pairwise(γ, domain₁, domain₂)

Evaluate variogram `γ` between all elements of `domain₁` and `domain₂`.
"""
function pairwise(γ::Variogram, domain₁, domain₂)
  u = first(domain₁)
  v = first(domain₂)
  m = nelements(domain₁)
  n = nelements(domain₂)
  R = result_type(γ, u, v)
  Γ = Matrix{R}(undef, m, n)
  pairwise!(Γ, γ, domain₁, domain₂)
end

function pairwise!(Γ, γ::Variogram, domain₁, domain₂)
  m = nelements(domain₁)
  n = nelements(domain₂)
  @inbounds for j in 1:n
    vj = domain₂[j]
    for i in 1:m
      ui = domain₁[i]
      Γ[i,j] = γ(ui, vj)
    end
  end
  Γ
end
