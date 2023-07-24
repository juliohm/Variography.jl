# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CressieEstimator()

Cressie's variogram estimator (equation 6.1 of Webster, R and Oliver, MA).

## References

* Webster, R and Oliver, MA. 2007. [Geostatistics for Environmental Scientists]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9780470517277)
"""
struct CressieEstimator <: VariogramEstimator end

function k(n)
  a, b, c = 0.457, 0.494, 0.045
  2 * n * n * (n * (n * a + b) + c)
end

formula(::CressieEstimator, z₁ᵢ, z₁ⱼ, z₂ᵢ, z₂ⱼ) = ((z₁ᵢ - z₁ⱼ) ⋅ (z₂ᵢ - z₂ⱼ))^(1 / 4)

normsum(::CressieEstimator, Σy, n) = Σy^4 / k(n)

combine(::CressieEstimator, yα, nα, yβ, nβ) = ((yα * k(nα))^(1 / 4) + (yβ * k(nβ))^(1 / 4))^4 / k(nα + nβ)
