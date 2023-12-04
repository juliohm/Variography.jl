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
  * estimator - variogram estimator (default to `:matheron` estimator)
  * algorithm - accumulation algorithm (default to `:ball`)

Available estimators:

  * `:matheron` - simple estimator based on squared differences
  * `:cressie`  - robust estimator based on 4th power of differences

Available algorithms:

  * `:full` - loop over all pairs of points in the data
  * `:ball` - loop over all points inside maximum lag ball

All implemented algorithms produce the exact same result.
The `:ball` algorithm is considerably faster when the
maximum lag is much smaller than the bounding box of
the domain of the data.

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
  data::AbstractGeoTable,
  var₁,
  var₂=var₁;
  nlags=20,
  maxlag=0.1diagonal(boundingbox(domain(data))),
  distance=Euclidean(),
  estimator=:matheron,
  algorithm=:ball
)

  # retrieve table and domain
  𝒯 = values(data)
  𝒟 = domain(data)

  # retrieve number of elements
  nelem = nelements(𝒟)

  # sanity checks
  @assert nelem > 1 "variogram requires at least 2 elements"
  @assert nlags > 0 "number of lags must be positive"
  @assert maxlag > 0 "maximum lag distance must be positive"
  @assert estimator ∈ (:matheron, :cressie) "invalid variogram estimator"
  @assert algorithm ∈ (:full, :ball) "invalid accumulation algorithm"

  # choose variogram estimator
  estim = estimator == :matheron ? MatheronEstimator() : CressieEstimator()

  # ball search with NearestNeighbors.jl requires AbstractFloat and MinkowskiMetric
  # https://github.com/KristofferC/NearestNeighbors.jl/issues/13
  isfloat = coordtype(𝒟) <: AbstractFloat
  isminkowski = distance isa MinkowskiMetric

  # warn users requesting :ball option with invalid parameters
  (algorithm == :ball && !isfloat) && @warn ":ball algorithm requires floating point coordinates, falling back to :full"
  (algorithm == :ball && !isminkowski) && @warn ":ball algorithm requires Minkowski metric, falling back to :full"

  # choose accumulation algorithm
  algo = if algorithm == :ball && isfloat && isminkowski
    BallSearchAccum(nlags, maxlag, distance)
  else
    FullSearchAccum(nlags, maxlag, distance)
  end

  # empirical variograms are defined on point sets
  𝒫 = PointSet([centroid(𝒟, i) for i in 1:nelem])
  𝒮 = georef(𝒯, 𝒫)

  # accumulate data with chosen algorithm
  abscissa, ordinate, counts = accumulate(𝒮, var₁, var₂, estim, algo)

  EmpiricalVariogram(abscissa, ordinate, counts, distance, estim)
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
assuming that both variograms have the same number of lags, distance
and estimator.
"""
function merge(γα::EmpiricalVariogram{V,D,E}, γβ::EmpiricalVariogram{V,D,E}) where {V,D,E}
  xα = γα.abscissa
  xβ = γβ.abscissa
  yα = γα.ordinate
  yβ = γβ.ordinate
  nα = γα.counts
  nβ = γβ.counts

  # copy distance and estimator
  d = γα.distance
  e = γα.estimator

  # merge function for estimator
  mergefun(yα, nα, yβ, nβ) = combine(e, yα, nα, yβ, nβ)

  # merge coordinates and bin counts
  n = nα + nβ
  x = @. (xα * nα + xβ * nβ) / n
  y = @. mergefun(yα, nα, yβ, nβ)

  # adjust empty bins
  x[n .== 0] .= xα[n .== 0]
  y[n .== 0] .= 0

  EmpiricalVariogram(x, y, n, d, e)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, γ::EmpiricalVariogram)
  ioctx = IOContext(io, :compact => true)
  print(ioctx, "EmpiricalVariogram(")
  print(ioctx, "abscissa: ")
  _printvec(ioctx, γ.abscissa, 1)
  print(ioctx, ", ordinate: ")
  _printvec(ioctx, γ.ordinate, 1)
  print(ioctx, ", distance: ", γ.distance)
  print(ioctx, ", estimator: ", γ.estimator)
  print(ioctx, ", npairs: ", sum(γ.counts))
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", γ::EmpiricalVariogram)
  ioctx = IOContext(io, :compact => true, :limit => true)
  println(ioctx, "EmpiricalVariogram")
  print(ioctx, "├─ abscissa: ")
  _printlnvec(ioctx, γ.abscissa, 3)
  print(ioctx, "├─ ordinate: ")
  _printlnvec(ioctx, γ.ordinate, 3)
  println(ioctx, "├─ distance: ", γ.distance)
  println(ioctx, "├─ estimator: ", γ.estimator)
  print(ioctx, "└─ npairs: ", sum(γ.counts))
end

function _printlnvec(io, vec, n)
  _printvec(io, vec, n)
  println(io)
end

function _printvec(io, vec, n)
  print(io, "[")
  if length(vec) > 2n
    k = n - 1
    join(io, vec[begin:(begin + k)], ", ")
    print(io, ", ..., ")
    join(io, vec[(end - k):end], ", ")
  else
    join(io, vec, ", ")
  end
  print(io, "]")
end
