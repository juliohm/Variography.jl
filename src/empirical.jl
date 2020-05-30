# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogram(sdata, var₁, var₂=var₁; [optional parameters])

Computes the empirical (a.k.a. experimental) omnidirectional
(cross-)variogram for variables `var₁` and `var₂` stored in
spatial data `sdata`.

    EmpiricalVariogram(partition, var₁, var₂=var₁; [optional parameters])

Alternatively, compute the (cross-)variogram on a `partition` of the data.

## Parameters

  * nlags    - number of lags (default to `20`)
  * maxlag   - maximum lag (default to half of maximum lag of data)
  * distance - custom distance function (default to `Euclidean` distance)
  * algo     - accumulation algorithm (default to `:ball`)

Available algorithms:

  * `:full` - loop over all pairs of points in the data
  * `:ball` - loop over all points inside maximum lag ball

All implemented algorithms produce the exact same result.
The `:ball` algorithm is considerably faster when the
maximum lag is much smaller than the bounding box of
the data.

See also: [`DirectionalVariogram`](@ref)

## References

* Chilès, JP and Delfiner, P. 2012. [Geostatistics: Modeling Spatial Uncertainty]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9781118136188)
"""
struct EmpiricalVariogram
  abscissa::Vector{Float64}
  ordinate::Vector{Float64}
  counts::Vector{Int}
end

function EmpiricalVariogram(sdata::AbstractData{T,N},
                            var₁::Symbol, var₂::Symbol=var₁;
                            nlags=20, maxlag=nothing,
                            distance=Euclidean(),
                            algo=:ball) where {N,T}
  # compute relevant parameters
  npts = npoints(sdata)
  hmax = isnothing(maxlag) ? 0.1diagonal(boundbox(sdata)) : maxlag

  # sanity checks
  @assert (var₁, var₂) ⊆ keys(variables(sdata)) "invalid variable names"
  @assert algo ∈ (:full, :ball) "invalid accumulation algorithm"
  @assert nlags > 0 "number of lags must be positive"
  @assert npts  > 1 "variogram requires at least 2 points"
  @assert hmax  > 0 "maximum lag distance must be positive"

  # ball search with NearestNeighbors.jl requires AbstractFloat
  # https://github.com/KristofferC/NearestNeighbors.jl/issues/13
  isfloat = coordtype(sdata) <: AbstractFloat

  # warn users requesting :ball option with non-floating point coordinates
  (algo == :ball && !isfloat) && @warn ":ball option requires floating point coordinates"

  # choose accumulation algorithm
  if isfloat && algo == :ball
    sums, counts = ball_search_accum(sdata, var₁, var₂, hmax, nlags, distance)
  else
    sums, counts = full_search_accum(sdata, var₁, var₂, hmax, nlags, distance)
  end

  # bin (or lag) size
  δh = hmax / nlags

  # variogram abscissa
  abscissa = range(δh/2, stop=hmax - δh/2, length=nlags)

  # variogram ordinate
  ordinate = @. (sums / counts) / 2
  ordinate[counts .== 0] .= NaN

  EmpiricalVariogram(abscissa, ordinate, counts)
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

# ------------------------
# ACCUMULATION ALGORITHMS
# ------------------------
function full_search_accum(sdata::AbstractData{T,N},
                           var₁::Symbol, var₂::Symbol,
                           hmax::T, nlags::Integer,
                           distance::Metric) where {N,T}
  # number of points to loop over
  npts = npoints(sdata)
  δh = hmax / nlags

  # lag sums and counts
  sums   = zeros(nlags)
  counts = zeros(Int, nlags)

  # preallocate memory for coordinates
  xi = MVector{N,T}(undef)
  xj = MVector{N,T}(undef)

  # loop over all pairs of points
  @inbounds for j in 1:npts
    coordinates!(xj, sdata, j)
    for i in j+1:npts
      coordinates!(xi, sdata, i)

      # evaluate spatial lag
      h = evaluate(distance, xi, xj)
      h > hmax && continue # early exit if out of range

      # evaluate (cross-)variance
      δ₁ = sdata[i,var₁] - sdata[j,var₁]
      δ₂ = sdata[i,var₂] - sdata[j,var₂]
      v = δ₁*δ₂

      # bin (or lag) where to accumulate result
      lag = ceil(Int, h / δh)

      if lag ≤ nlags && !ismissing(v) && !isnan(v)
        sums[lag] += v
        counts[lag] += 1
      end
    end
  end

  sums, counts
end

function ball_search_accum(sdata::AbstractData{T,N},
                           var₁::Symbol, var₂::Symbol,
                           hmax::T, nlags::Integer,
                           distance::Metric) where {N,T}
  # retrieve relevant parameters
  npts = npoints(sdata)
  δh = hmax / nlags

  # lag sums and counts
  sums   = zeros(nlags)
  counts = zeros(Int, nlags)

  # preallocate memory for coordinates
  xi = MVector{N,T}(undef)
  xj = MVector{N,T}(undef)

  # fast ball search
  ball = BallNeighborhood{N}(hmax, distance)
  searcher = NeighborhoodSearcher(sdata, ball)

  # loop over points inside norm ball
  @inbounds for j in 1:npts
    coordinates!(xj, sdata, j)
    for i in search(xj, searcher)
      i ≤ j && continue # avoid double counting
      coordinates!(xi, sdata, i)

      # evaluate spatial lag
      h = evaluate(distance, xi, xj)

      # evaluate (cross-)variance
      δ₁ = sdata[i,var₁] - sdata[j,var₁]
      δ₂ = sdata[i,var₂] - sdata[j,var₂]
      v = δ₁*δ₂

      # bin (or lag) where to accumulate result
      lag = ceil(Int, h / δh)

      if lag ≤ nlags && !ismissing(v) && !isnan(v)
        sums[lag] += v
        counts[lag] += 1
      end
    end
  end

  sums, counts
end
