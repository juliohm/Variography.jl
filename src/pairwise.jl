# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    pairwise(γ, domain)

Evaluate variogram `γ` between all elements in the `domain`.
"""
function pairwise(γ::Variogram, domain)
  u = first(domain)
  n = length(domain)
  R = result_type(γ, u, u)
  Γ = Matrix{R}(undef, n, n)
  pairwise!(Γ, γ, domain)
end

function pairwise!(Γ, γ::Variogram, domain)
  n = length(domain)
  @inbounds for j in 1:n
    vⱼ = domain[j]
    sⱼ = _sample(γ, vⱼ)
    for i in (j + 1):n
      vᵢ = domain[i]
      sᵢ = _sample(γ, vᵢ)
      Γ[i, j] = mean(γ(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ)
    end
    Γ[j, j] = mean(γ(pⱼ, pⱼ) for pⱼ in sⱼ, pⱼ in sⱼ)
    for i in 1:(j - 1)
      Γ[i, j] = Γ[j, i] # leverage the symmetry
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
  m = length(domain₁)
  n = length(domain₂)
  R = result_type(γ, u, v)
  Γ = Matrix{R}(undef, m, n)
  pairwise!(Γ, γ, domain₁, domain₂)
end

function pairwise!(Γ, γ::Variogram, domain₁, domain₂)
  m = length(domain₁)
  n = length(domain₂)
  @inbounds for j in 1:n
    vⱼ = domain₂[j]
    sⱼ = _sample(γ, vⱼ)
    for i in 1:m
      vᵢ = domain₁[i]
      sᵢ = _sample(γ, vᵢ)
      Γ[i, j] = mean(γ(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ)
    end
  end
  Γ
end
