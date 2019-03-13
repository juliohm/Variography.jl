# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogram(X, z₁, z₂=z₁; [optional parameters])

Computes the empirical (a.k.a. experimental) omnidirectional
(cross-)variogram from data locations `X` and values `z₁` and `z₂`.

    EmpiricalVariogram(spatialdata, var₁, var₂=var₁; [optional parameters])

Alternatively, compute the (cross-)variogram for the variables
`var₁` and `var₂` stored in a `spatialdata` object.

    EmpiricalVariogram(partition, var₁, var₂=var₁; [optional parameters])

Alternatively, compute the (cross-)variogram on a `partition` of the data.

## Parameters

  * nlags - number of lags (default to 20)
  * maxlag - maximum lag (default to half of maximum lag of data)
  * distance - custom distance function (default to Euclidean distance)

See also: [`DirectionalVariogram`](@ref)
"""
struct EmpiricalVariogram
  abscissa::Vector{Float64}
  ordinate::Vector{Float64}
  counts::Vector{Int}

  function EmpiricalVariogram(X, z₁, z₂, nlags, maxlag, distance)
    # sanity checks
    @assert nlags > 0 "number of lags must be positive"
    if maxlag ≠ nothing
      @assert maxlag > 0 "maximum lag distance must be positive"
    end

    # number of point pairs
    npoints = size(X, 2)
    npairs = (npoints * (npoints-1)) ÷ 2

    # compute pairwise distance
    lags  = Vector{Float64}(undef, npairs)
    zdiff = Vector{Float64}(undef, npairs)
    idx = 1
    for j=1:npoints
      xj = view(X,:,j)
      for i=j+1:npoints
        xi = view(X,:,i)
        @inbounds δx = evaluate(distance, xi, xj)
        @inbounds δz = (z₁[i] - z₁[j])*(z₂[i] - z₂[j])
        @inbounds lags[idx] = δx
        @inbounds zdiff[idx] = ismissing(δz) ? NaN : δz
        idx += 1
      end
    end

    # default maximum lag
    maxlag == nothing && (maxlag = maximum(lags) / 2)

    # corresponding bin size
    binsize = maxlag / nlags

    # variogram abscissa
    abscissa = range(binsize/2, stop=maxlag - binsize/2, length=nlags)

    # handle missing/invalid values
    valid = @. !isnan(zdiff)

    if any(valid)
      lags  = lags[valid]
      zdiff = zdiff[valid]

      # find bin for the pair
      binidx  = ceil.(Int, lags / binsize)

      # discard lags greater than maximum lag
      zdiff  = zdiff[binidx .≤ nlags]
      binidx = binidx[binidx .≤ nlags]

      # place squared differences at the bins
      bins = [zdiff[binidx .== i] for i=1:nlags]

      # variogram ordinate and count
      ordinate = [length(bin) > 0 ? sum(bin)/2length(bin) : NaN for bin in bins]
      counts   = length.(bins)

      new(abscissa, ordinate, counts)
    else
      new(abscissa, fill(NaN, nlags), fill(0, nlags))
    end
  end
end

EmpiricalVariogram(X, z₁, z₂=z₁; nlags=20, maxlag=nothing, distance=Euclidean()) =
  EmpiricalVariogram(X, z₁, z₂, nlags, maxlag, distance)

function EmpiricalVariogram(spatialdata::S, var₁::Symbol, var₂::Symbol=var₁;
                            kwargs...) where {S<:AbstractSpatialData}
  npts = npoints(spatialdata)

  X = coordinates(spatialdata)
  z₁ = [value(spatialdata, i, var₁) for i in 1:npts]
  z₂ = var₁ ≠ var₂ ? [value(spatialdata, i, var₂) for i in 1:npts] : z₁

  EmpiricalVariogram(X, z₁, z₂; kwargs...)
end

function EmpiricalVariogram(partition::SpatialPartition,
                            var₁::Symbol, var₂::Symbol=var₁; kwargs...)
  # retain spatial data with at least 2 points
  filtered = Iterators.filter(d -> npoints(d) > 1, partition)

  @assert !isempty(filtered) "invalid partition of spatial data"

  spatialdata, _ = iterate(filtered)
  γ = EmpiricalVariogram(spatialdata, var₁, var₂; kwargs...)
  for spatialdata in Iterators.drop(filtered, 1)
    γiter = EmpiricalVariogram(spatialdata, var₁, var₂; kwargs...)
    merge!(γ, γiter)
  end

  γ
end

"""
    DirectionalVariogram(spatialdata, direction, var₁, var₂=var₁; [optional parameters])

Computes the empirical (cross-)variogram for the variables `var₁` and `var₂` stored in
`spatialdata` along a given `direction`.

Optional parameters include the parameters for `EmpiricalVariogram` and the parameters
for `DirectionPartitioner`.

### Notes

A `DirectionalVariogram` is just a function that first partitions the `spatialdata`
using a `DirectionPartitioner` and then passes the result to the corresponding
`EmpiricalVariogram` constructor.

See also: [`EmpiricalVariogram`](@ref), [`DirectionPartitioner`](@ref)
"""
function DirectionalVariogram(spatialdata::S, direction::NTuple,
                              var₁::Symbol, var₂::Symbol=var₁;
                              tol=1e-6, kwargs...) where {S<:AbstractSpatialData}
  partitioner = DirectionPartitioner(direction; tol=tol)
  EmpiricalVariogram(partition(spatialdata, partitioner), var₁, var₂; kwargs...)
end

"""
    merge!(γₐ, γᵦ)

Merge the empirical variogram `γₐ` with the empirical variogram `γᵦ`
assuming that both have the same abscissa.
"""
function merge!(γₐ::EmpiricalVariogram, γᵦ::EmpiricalVariogram)
  yₐ = γₐ.ordinate
  yᵦ = γᵦ.ordinate
  nₐ = γₐ.counts
  nᵦ = γᵦ.counts

  for i in eachindex(yₐ)
    if nᵦ[i] > 0
      if nₐ[i] > 0
        yₐ[i] = (yₐ[i]*nₐ[i] + yᵦ[i]*nᵦ[i]) / (nₐ[i] + nᵦ[i])
      else
        yₐ[i] = yᵦ[i]
      end
    end
  end

  @. nₐ = nₐ + nᵦ

  nothing
end

"""
    values(γ)

Returns the center of the bins, the mean squared differences divided by 2
and the number of squared differences at the bins for a given empirical
variogram `γ`.

## Examples

Plotting empirical variogram manually:

```julia
julia> x, y, n = values(γemp)
julia> plot(x, y, label="variogram")
julia> bar!(x, n, label="histogram")
```
"""
Base.values(γ::EmpiricalVariogram) = γ.abscissa, γ.ordinate, γ.counts
