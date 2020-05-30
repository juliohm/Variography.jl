# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogram(partition, var₁, var₂=var₁; [parameters])

Compute the empirical (cross-)variogram of the spatial `partition` for
variables `var₁` and `var₂`.

## References

* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with partition variograms]
  (https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
function EmpiricalVariogram(partition::SpatialPartition,
                            var₁::Symbol, var₂::Symbol=var₁; kwargs...)
  # retain spatial data with at least 2 points
  filtered = Iterators.filter(d -> npoints(d) > 1, partition)

  @assert !isempty(filtered) "invalid partition of spatial data"

  sdata, _ = iterate(filtered)
  γ = EmpiricalVariogram(sdata, var₁, var₂; kwargs...)
  for sdata in Iterators.drop(filtered, 1)
    γiter = EmpiricalVariogram(sdata, var₁, var₂; kwargs...)
    merge!(γ, γiter)
  end

  γ
end

"""
    DirectionalVariogram(sdata, direction, var₁, var₂=var₁; dtol=1e-6, [parameters])

Computes the empirical (cross-)variogram for the variables `var₁` and `var₂` stored in
spatial data `sdata` along a given `direction`.

Optional parameters include the parameters for [`EmpiricalVariogram`](@ref) and the
parameters for [`DirectionPartitioner`](@ref).
"""
function DirectionalVariogram(sdata::S, direction::NTuple,
                              var₁::Symbol, var₂::Symbol=var₁;
                              dtol=1e-6, kwargs...) where {S<:AbstractData}
  p = partition(sdata, DirectionPartitioner(direction; tol=dtol))
  EmpiricalVariogram(p, var₁, var₂; kwargs...)
end
