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

## Parameters

  * nlags - number of lags (default to 20)
  * maxlag - maximum lag (default to maximum lag of data)
  * distance - custom distance function (default to Euclidean distance)
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
    R = result_type(distance, view(X,:,1), view(X,:,1))

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
    valid = .!(ismissing.(zdiff) .| isnan.(zdiff))
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
  npts = npoints(spatialdata)

  X = hcat([coordinates(spatialdata, i) for i in 1:npts]...)
  z₁ = [value(spatialdata, i, var₁) for i in 1:npts]
  z₂ = var₁ ≠ var₂ ? [value(spatialdata, i, var₂) for i in 1:npts] : z₁

  EmpiricalVariogram(X, z₁, z₂; kwargs...)
end

"""
    values(γemp)

Returns the center of the bins, the mean squared differences divided by 2
and the number of squared differences at the bins for a given empirical
variogram `γemp`.

## Examples

Plotting empirical variogram manually:

```julia
julia> x, y, n = values(γemp)
julia> plot(x, y, label="variogram")
julia> bar!(x, n, label="histogram")
```
"""
Base.values(γ::EmpiricalVariogram) = γ.abscissa, γ.ordinate, γ.counts
