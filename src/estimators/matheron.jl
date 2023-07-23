# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MatheronEstimator()

Matheron's variogram estimator.
"""
struct MatheronEstimator <: VariogramEstimator end

formula(::MatheronEstimator, z₁ᵢ, z₁ⱼ, z₂ᵢ, z₂ⱼ) = (z₁ᵢ - z₁ⱼ) ⋅ (z₂ᵢ - z₂ⱼ)