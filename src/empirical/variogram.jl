# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogram(data, var₁, var₂=var₁; [parameters])

Computes the empirical (a.k.a. experimental) omnidirectional
(cross-)variogram for variables `var₁` and `var₂` stored in
geospatial `data`.

## Parameters

  * nlags     - number of lags (default to `20`)
  * maxlag    - maximum lag (default to 1/10 of maximum lag of data)
  * distance  - custom distance function (default to `Euclidean` distance)
  * estimator - variogram estimator (default to `Matheron` estimator)
  * algo      - accumulation algorithm (default to `:ball`)

Available algorithms:

  * `:full` - loop over all pairs of points in the data
  * `:ball` - loop over all points inside maximum lag ball

All implemented algorithms produce the exact same result.
The `:ball` algorithm is considerably faster when the
maximum lag is much smaller than the bounding box of
the data.

See also: [`DirectionalVariogram`](@ref), [`PlanarVariogram`](@ref),
[`EmpiricalVarioplane`](@ref).

## References

* Chilès, JP and Delfiner, P. 2012. [Geostatistics: Modeling Spatial Uncertainty]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9781118136188)
* Webster, R and Oliver, MA. 2007. [Geostatistics for Environmental Scientists]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9780470517277)
* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with partition variograms]
  (https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
struct EmpiricalVariogram{V,D,E}
  abscissa::Vector{Float64}
  ordinate::Vector{V}
  counts::Vector{Int}
  distance::D
  estimator::E
end

function EmpiricalVariogram(
  data,
  var₁::Symbol,
  var₂::Symbol=var₁;
  nlags=20,
  maxlag=0.1diagonal(boundingbox(data)),
  distance=Euclidean(),
  estimator=Matheron(),
  algo=:ball
)

  # retrieve table and domain
  𝒯 = values(data)
  𝒟 = domain(data)

  # retrieve the column names of data values
  vars = Tables.columnnames(𝒯)

  # sanity checks
  @assert nelements(𝒟) > 1 "variogram requires at least 2 elements"
  @assert (var₁, var₂) ⊆ vars "invalid variable names"
  @assert algo ∈ (:full, :ball) "invalid accumulation algorithm"
  @assert nlags > 0 "number of lags must be positive"
  @assert maxlag > 0 "maximum lag distance must be positive"

  # ball search with NearestNeighbors.jl requires AbstractFloat and MinkowskiMetric
  # https://github.com/KristofferC/NearestNeighbors.jl/issues/13
  isfloat = coordtype(𝒟) <: AbstractFloat
  isminkowski = distance isa MinkowskiMetric

  # warn users requesting :ball option with invalid parameters
  (algo == :ball && !isfloat) && @warn ":ball algorithm requires floating point coordinates"
  (algo == :ball && !isminkowski) && @warn ":ball algorithm requires Minkowski metric"

  # empirical variograms are defined on point sets
  𝒫 = PointSet([centroid(𝒟, i) for i in 1:nelements(𝒟)])
  𝒮 = georef(𝒯, 𝒫)

  # choose accumulation algorithm
  if algo == :ball && isfloat && isminkowski
    xsums, ysums, counts = ball_search_accum(𝒮, var₁, var₂, maxlag, nlags, distance, estimator)
  else
    xsums, ysums, counts = full_search_accum(𝒮, var₁, var₂, maxlag, nlags, distance, estimator)
  end

  # bin (or lag) size
  δh = maxlag / nlags
  lags = range(δh / 2, stop=maxlag - δh / 2, length=nlags)

  # variogram abscissa
  abscissa = @. xsums / counts
  abscissa[counts .== 0] .= lags[counts .== 0]

  # variogram ordinate
  ordinate = @. (ysums / counts) / 2
  ordinate[counts .== 0] .= zero(eltype(ordinate))

  EmpiricalVariogram(abscissa, ordinate, counts, distance, estimator)
end

"""
    values(γ)

Returns the abscissa, the ordinate, and the bin counts
of the empirical variogram `γ`.
"""
Base.values(γ::EmpiricalVariogram) = γ.abscissa, γ.ordinate, γ.counts

"""
    distance(γ)

Return the distance used to compute the empirical variogram `γ`.
"""
distance(γ::EmpiricalVariogram) = γ.distance

"""
    estimator(γ)

Return the estimator used to compute the empirical variogram `γ`.
"""
estimator(γ::EmpiricalVariogram) = γ.estimator

"""
    merge(γα, γβ)

Merge the empirical variogram `γα` with the empirical variogram `γβ`
assuming that both variograms have the same number of lags.
"""
function merge(γα::EmpiricalVariogram{V,D,E}, γβ::EmpiricalVariogram{V,D,E}) where {V,D,E}
  xα = γα.abscissa
  xβ = γβ.abscissa
  yα = γα.ordinate
  yβ = γβ.ordinate
  nα = γα.counts
  nβ = γβ.counts

  d = γα.distance
  e = γα.estimator

  n = nα + nβ
  x = @. (xα * nα + xβ * nβ) / n
  y = @. (yα * nα + yβ * nβ) / n

  # adjust empty bins
  x[n .== 0] .= xα[n .== 0]
  y[n .== 0] .= 0

  EmpiricalVariogram(x, y, n, d, e)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, ::EmpiricalVariogram)
  print(io, "EmpiricalVariogram")
end

function Base.show(io::IO, ::MIME"text/plain", γ::EmpiricalVariogram)
  println(io, γ)
  println(io, "  abscissa: ", extrema(γ.abscissa))
  println(io, "  ordinate: ", extrema(γ.ordinate))
  print(io, "  N° pairs: ", sum(γ.counts))
end
