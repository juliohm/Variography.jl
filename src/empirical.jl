# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

# KD-tree search is only valid for Minkowski metric
const MinkowskiMetric = Union{Euclidean,Chebyshev,Cityblock,Minkowski,
                              WeightedEuclidean,WeightedCityblock,
                              WeightedMinkowski}

# default maximum lag to be used in various methods
default_maxlag(sdata) = 0.1diagonal(boundbox(sdata))

"""
    EmpiricalVariogram(sdata, var₁, var₂=var₁; [optional parameters])

Computes the empirical (a.k.a. experimental) omnidirectional
(cross-)variogram for variables `var₁` and `var₂` stored in
spatial data `sdata`.

    EmpiricalVariogram(partition, var₁, var₂=var₁; [optional parameters])

Alternatively, compute the (cross-)variogram on a `partition` of the data
as described in Hoffimann & Zadrozny 2019.

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
* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with partition variograms]
  (https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
struct EmpiricalVariogram{D}
  abscissa::Vector{Float64}
  ordinate::Vector{Float64}
  counts::Vector{Int}
  distance::D
end

EmpiricalVariogram(abscissa, ordinate, counts, distance) =
  EmpiricalVariogram{typeof(distance)}(abscissa, ordinate, counts, distance)

function EmpiricalVariogram(sdata, var₁::Symbol, var₂::Symbol=var₁;
                            nlags=20, maxlag=default_maxlag(sdata),
                            distance=Euclidean(), algo=:ball)
  # relevant parameters
  N = ncoords(sdata)
  T = coordtype(sdata)
  npts = nelms(sdata)
  hmax = maxlag

  # sanity checks
  @assert (var₁, var₂) ⊆ name.(variables(sdata)) "invalid variable names"
  @assert algo ∈ (:full, :ball) "invalid accumulation algorithm"
  @assert nlags > 0 "number of lags must be positive"
  @assert npts  > 1 "variogram requires at least 2 points"
  @assert hmax  > 0 "maximum lag distance must be positive"

  # ball search with NearestNeighbors.jl requires AbstractFloat and MinkowskiMetric
  # https://github.com/KristofferC/NearestNeighbors.jl/issues/13
  isfloat     = coordtype(sdata) <: AbstractFloat
  isminkowski = distance isa MinkowskiMetric

  # warn users requesting :ball option with invalid parameters
  (algo == :ball && !isfloat) && @warn ":ball algorithm requires floating point coordinates"
  (algo == :ball && !isminkowski) && @warn ":ball algorithm requires Minkowski metric"

  # choose accumulation algorithm
  if algo == :ball && isfloat && isminkowski
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
  ordinate[counts .== 0] .= 0

  EmpiricalVariogram(abscissa, ordinate, counts, distance)
end

"""
    values(γ)

Returns the abscissa, the ordinate, and the bin counts
of the empirical variogram `γ`.
"""
Base.values(γ::EmpiricalVariogram) = γ.abscissa, γ.ordinate, γ.counts

"""
    distance(γ)

Retun the distance used to compute the empirical variogram `γ`.
"""
distance(γ::EmpiricalVariogram) = γ.distance

"""
    merge(γα, γβ)

Merge the empirical variogram `γα` with the empirical variogram `γβ`
assuming that both variograms have the same abscissa.
"""
function merge(γα::EmpiricalVariogram{D}, γβ::EmpiricalVariogram{D}) where {D}
  yα = γα.ordinate
  yβ = γβ.ordinate
  nα = γα.counts
  nβ = γβ.counts

  n = nα + nβ
  x = γα.abscissa
  y = @. (yα*nα + yβ*nβ) / n
  y[n .== 0] .= 0

  d = γα.distance

  EmpiricalVariogram(x, y, n, d)
end

# ------------
# IO methods
# ------------
function Base.show(io::IO, γ::EmpiricalVariogram)
  print(io, "EmpiricalVariogram")
end

function Base.show(io::IO, ::MIME"text/plain", γ::EmpiricalVariogram)
  println(io, γ)
  println(io, "  abscissa: ", extrema(γ.abscissa))
  println(io, "  ordinate: ", extrema(γ.ordinate))
  println(io, "  N° pairs: ", sum(γ.counts))
end

# ------------------------
# ACCUMULATION ALGORITHMS
# ------------------------
function full_search_accum(sdata, var₁, var₂, hmax, nlags, distance)
  # retrieve relevant parameters
  N = ncoords(sdata)
  T = coordtype(sdata)
  npts = nelms(sdata)
  δh = hmax / nlags

  # lag sums and counts
  sums   = zeros(nlags)
  counts = zeros(Int, nlags)

  # preallocate memory for coordinates
  xi = MVector{N,T}(undef)
  xj = MVector{N,T}(undef)

  # collect vectors for variables
  Z₁, Z₂ = sdata[var₁], sdata[var₂]

  # loop over all pairs of points
  @inbounds for j in 1:npts
    coordinates!(xj, sdata, j)
    for i in j+1:npts
      coordinates!(xi, sdata, i)

      # evaluate spatial lag
      h = evaluate(distance, xi, xj)
      h > hmax && continue # early exit if out of range

      # evaluate (cross-)variance
      v = (Z₁[i] - Z₁[j]) * (Z₂[i] - Z₂[j])

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

function ball_search_accum(sdata, var₁, var₂, hmax, nlags, distance)
  # retrieve relevant parameters
  N = ncoords(sdata)
  T = coordtype(sdata)
  npts = nelms(sdata)
  δh = hmax / nlags

  # lag sums and counts
  sums   = zeros(nlags)
  counts = zeros(Int, nlags)

  # preallocate memory for coordinates
  xi = MVector{N,T}(undef)
  xj = MVector{N,T}(undef)

  # collect vectors for variables
  Z₁, Z₂ = sdata[var₁], sdata[var₂]

  # fast ball search
  ball = BallNeighborhood(hmax, distance)
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
      v = (Z₁[i] - Z₁[j]) * (Z₂[i] - Z₂[j])

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