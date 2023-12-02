# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# models that can be fitted currently
fittable() = filter(isstationary, setdiff(subtypes(Variogram), (NuggetEffect, NestedVariogram)))

"""
    VariogramFitAlgo

An algorithm for fitting theoretical variograms.
"""
abstract type VariogramFitAlgo end

"""
    WeightedLeastSquares()
    WeightedLeastSquares(w)

Fit theoretical variogram using weighted least squares with weighting
function `w` (e.g. h -> 1/h). If no weighting function is provided,
bin counts of empirical variogram are normalized and used as weights.
"""
struct WeightedLeastSquares <: VariogramFitAlgo
  weightfun::Union{Function,Nothing}
end

WeightedLeastSquares() = WeightedLeastSquares(nothing)

"""
    fit(V, g, algo=WeightedLeastSquares())

Fit theoretical variogram type `V` to empirical variogram `g`
using algorithm `algo`.

## Examples

```julia
julia> fit(SphericalVariogram, g)
julia> fit(ExponentialVariogram, g)
julia> fit(GaussianVariogram, g, WeightedLeastSquares())
```
"""
fit(V::Type{<:Variogram}, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares()) =
  fit_impl(V, g, algo) |> first

"""
    fit(Vs, g, algo=WeightedLeastSquares())

Fit theoretical variogram types `Vs` to empirical variogram `g`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit([SphericalVariogram, ExponentialVariogram], g)
```
"""
function fit(Vs, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares())
  # fit each variogram type
  res = [fit_impl(V, g, algo) for V in Vs]
  γs, ϵs = first.(res), last.(res)

  # return best candidate
  γs[argmin(ϵs)]
end

"""
    fit(Variogram, g, algo=WeightedLeastSquares())

Fit all "fittable" subtypes of `Variogram` to empirical variogram `g`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit(Variogram, g)
julia> fit(Variogram, g, WeightedLeastSquares())
```

See also `Variography.fittable()`.
"""
fit(::Type{Variogram}, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares()) = fit(fittable(), g, algo)

"""
    fit(V, g, weightfun)
    fit(Variogram, g, weightfun)

Convenience method that forwards the weighting function `weightfun`
to the `WeightedLeastSquares` algorithm.

## Examples

```julia
fit(SphericalVariogram, g, h -> exp(-h))
fit(Variogram, g, h -> exp(-h/100))
```
"""
fit(V, g::EmpiricalVariogram, weightfun::Function) = fit(V, g, WeightedLeastSquares(weightfun))

# ---------------
# IMPLEMENTATION
# ---------------

function fit_impl(V::Type{<:Variogram}, g::EmpiricalVariogram, algo::WeightedLeastSquares)
  # values of empirical variogram
  x, y, n = values(g)

  # custom ball of given radius
  ball(r) = MetricBall(r, distance(g))

  # discard invalid bins
  x = x[n .> 0]
  y = y[n .> 0]
  n = n[n .> 0]

  # strip units if necessary
  𝓊 = unit(first(y))
  y = ustrip.(y)

  # auxiliary variables
  xmax, ymax = maximum(x), maximum(y)

  # evaluate weights
  f = algo.weightfun
  w = isnothing(f) ? n / sum(n) : map(f, x)

  # objective function
  function J(θ)
    γ = V(ball(θ[1]), sill=θ[2] + θ[3], nugget=θ[3])
    sum(w[i] * (γ(x[i]) - y[i])^2 for i in eachindex(x))
  end

  # initial guess
  θₒ = [xmax / 3, 0.95 * ymax, 1e-6]

  # box constraints
  l = [0.0, 0.0, 0.0]
  u = [xmax, ymax, ymax]

  # solve optimization problem
  sol = Optim.optimize(J, l, u, θₒ)
  ϵ = Optim.minimum(sol)
  θ = Optim.minimizer(sol)

  # optimal variogram (with units)
  γ = V(ball(θ[1]), sill=(θ[2] + θ[3]) * 𝓊, nugget=θ[3] * 𝓊)

  γ, ϵ
end
