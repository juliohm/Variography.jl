# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MatheronEstimator()

Matheron's variogram estimator.
"""
struct MatheronEstimator <: VariogramEstimator end

result_type(::MatheronEstimator, z₁, z₂) = typeof((z₁[1] - z₂[1]) ⋅ (z₁[1] - z₂[1]))

formula(::MatheronEstimator, z₁ᵢ, z₁ⱼ, z₂ᵢ, z₂ⱼ) = (z₁ᵢ - z₁ⱼ) ⋅ (z₂ᵢ - z₂ⱼ)