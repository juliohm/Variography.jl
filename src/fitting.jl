# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

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
    fit(V, γ, [algo])

Fit theoretical variogram type `V` to empirical variogram `γ`
using algorithm `algo`. Default algorithm is `WeightedLeastSquares`.

## Examples

```julia
julia> fit(SphericalVariogram, γ)
```
"""
function fit(V::Type{<:Variogram}, γ::EmpiricalVariogram,
             algo::VariogramFitAlgo=WeightedLeastSquares())
  # dispatch appropriate implementation
  vario, err = fit_impl(V, γ, algo)

  vario
end

"""
    fit(Variogram, γ, [algo])

Fit all stationary variogram types to empirical variogram `γ`,
which are subtypes of `Variogram`, and return the one with
minimum error as defined by the algorithm `algo`.

## Examples

```julia
julia> fit(Variogram, γ)
```
"""
function fit(::Type{Variogram}, γ::EmpiricalVariogram,
             algo::VariogramFitAlgo=WeightedLeastSquares())
  # list of variogram types to try
  Vs = filter(isstationary, setdiff(subtypes(Variogram), [CompositeVariogram]))

  # fit each variogram type
  res = [fit_impl(V, γ, algo) for V in Vs]
  γs, es = first.(res), last.(res)

  # return best candidate
  γs[argmin(es)]
end

function fit_impl(V::Type{<:Variogram}, γ::EmpiricalVariogram,
                  algo::WeightedLeastSquares)
  # retrieve empirical points
  x, y, n = values(γ)

  # discard invalid bins
  x = x[n .> 0]
  y = y[n .> 0]
  n = n[n .> 0]

  # evaluate weights
  f = algo.weightfun
  w = f ≠ nothing ? map(f, x) : n / sum(n)

  # objective function
  J(p) = begin
    g = V(range=p[1], sill=p[2], nugget=p[3])
    sum(w[i]*(g(x[i]) - y[i])^2 for i in eachindex(x))
  end

  # auxiliary variables
  xmax, ymax = maximum(x), maximum(y)

  # initial guess
  pₒ = [xmax/3, .95*ymax, 1e-6]

  # box constraints
  l  = [0., 0., 0.]
  u  = [xmax, ymax, ymax]

  # solve optimization problem
  sol = Optim.optimize(J, l, u, pₒ)
  err = Optim.minimum(sol)
  p   = Optim.minimizer(sol)

  # best variogram
  vario = V(range=p[1], sill=p[2], nugget=p[3])

  vario, err
end
