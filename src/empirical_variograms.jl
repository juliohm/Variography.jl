# ------------------------------------------------------------------
# Copyright (c) 2017, Júlio Hoffimann Mendes <juliohm@stanford.edu>
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
  * maxlag - maximum lag (default to maximum lag of data)
  * distance - custom distance function (default to Euclidean distance)

See also: [`DirectionalVariogram`](@ref)
"""
struct EmpiricalVariogram{T<:Real,V,D<:Metric}
  abscissa::Vector{Float64}
  ordinate::Vector{Float64}
  counts::Vector{Int}

  function EmpiricalVariogram{T,V,D}(X, z₁, z₂,
                                     nlags, maxlag,
                                     distance) where {T<:Real,V,D<:Metric}
    # sanity checks
    @assert nlags > 0 "number of lags must be positive"
    if maxlag ≠ nothing
      @assert maxlag > 0 "maximum lag distance must be positive"
    end

    # number of point pairs
    npoints = size(X, 2)
    npairs = (npoints * (npoints-1)) ÷ 2

    # result type of distance between coordinates
    R = Distances.result_type(distance, view(X,:,1), view(X,:,1))

    # compute pairwise distance
    lags  = Vector{R}(undef, npairs)
    zdiff = Vector{V}(undef, npairs)
    idx = 1
    for j=1:npoints
      xj = view(X,:,j)
      for i=j+1:npoints
        xi = view(X,:,i)
        @inbounds lags[idx] = evaluate(distance, xi, xj)
        @inbounds zdiff[idx] = (z₁[i] - z₁[j])*(z₂[i] - z₂[j])
        idx += 1
      end
    end

    # handle missing/invalid values
    valid = @. !(ismissing(zdiff) | isnan(zdiff))
    lags  = lags[valid]
    zdiff = zdiff[valid]

    # default maximum lag
    maxlag == nothing && (maxlag = maximum(lags))

    # find bin for the pair
    binsize = maxlag / nlags
    binidx  = ceil.(Int, lags / binsize)

    # discard lags greater than maximum lag
    zdiff  = zdiff[binidx .≤ nlags]
    binidx = binidx[binidx .≤ nlags]

    # place squared differences at the bins
    bins = [zdiff[binidx .== i] for i=1:nlags]

    # variogram abscissa, ordinate, and count
    abscissa = range(binsize/2, stop=maxlag - binsize/2, length=nlags)
    ordinate = [length(bin) > 0 ? sum(bin)/2length(bin) : NaN for bin in bins]
    counts   = length.(bins)

    new(abscissa, ordinate, counts)
  end
end

EmpiricalVariogram(X, z₁, z₂=z₁; nlags=20, maxlag=nothing, distance=Euclidean()) =
  EmpiricalVariogram{eltype(X),eltype(z₁),typeof(distance)}(X, z₁, z₂, nlags, maxlag, distance)

function EmpiricalVariogram(spatialdata::S, var₁::Symbol, var₂::Symbol=var₁;
                            kwargs...) where {S<:AbstractSpatialData}
  T = coordtype(spatialdata)
  N = ndims(spatialdata)
  npts = npoints(spatialdata)

  X = Matrix{T}(undef, N, npts)
  for i in 1:npts
    coordinates!(view(X,:,i), spatialdata, i)
  end
  z₁ = [value(spatialdata, i, var₁) for i in 1:npts]
  z₂ = var₁ ≠ var₂ ? [value(spatialdata, i, var₂) for i in 1:npts] : z₁

  EmpiricalVariogram(X, z₁, z₂; kwargs...)
end

function EmpiricalVariogram(partition::P, var₁::Symbol, var₂::Symbol=var₁;
                            kwargs...) where {P<:AbstractPartition}
  spatialdata, state = iterate(partition)
  γ = EmpiricalVariogram(spatialdata, var₁, var₂; kwargs...)
  while state < length(partition)
    spatialdata, state = iterate(partition, state)
    γiter = EmpiricalVariogram(spatialdata, var₁, var₂; kwargs...)
    update!(γ, γiter)
  end

  γ
end

"""
    DirectionalVariogram(spatialdata, direction, var₁, var₂=var₁; [optional parameters])

Computes the empirical (cross-)variogram for the variables `var₁` and `var₂` stored in
`spatialdata` along a given `direction`.

### Notes

A `DirectionalVariogram` is just a function that first partitions the `spatialdata`
using a `DirectionalPartition` and then passes the result to the corresponding
`EmpiricalVariogram` constructor.

See also: [`EmpiricalVariogram`](@ref)
"""
function DirectionalVariogram(spatialdata::S, direction::NTuple,
                              var₁::Symbol, var₂::Symbol=var₁;
                              kwargs...) where {S<:AbstractSpatialData}
  EmpiricalVariogram(DirectionalPartition(spatialdata, direction), var₁, var₂; kwargs...)
end

"""
    update!(γₐ, γᵦ)

Update the empirical variogram `γₐ` with the empirical variogram `γᵦ`
assuming that both have the same abscissa.
"""
function update!(γₐ::EmpiricalVariogram{T,V,D},
                 γᵦ::EmpiricalVariogram{T,V,D}) where {T<:Real,V,D<:Metric}
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
