# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# models that can be fitted currently
const FITTABLE = filter(isstationary, setdiff(subtypes(Variogram), (NuggetEffect, NestedVariogram)))

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
    fit(V, γ, [weightfun])

Fit theoretical variogram type `V` to empirical variogram `γ`
using algorithm `algo`. Default algorithm is `WeightedLeastSquares`.

Alternatively pass the weighting function `weightfun` directly
to the fitting procedure.
"""
function fit(V::Type{<:Variogram}, γ::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares())
  # dispatch appropriate implementation
  vario, err = fit_impl(V, γ, algo)

  vario
end

"""
    fit(Variogram, γ, [algo])
    fit(Variogram, γ, [weightfun])

Fit all subtypes of `Variogram` to empirical variogram `γ` and
return the one with minimum error as defined by the algorithm `algo`.

Alternatively pass the weighting function `weightfun` directly
to the fitting procedure.
"""
function fit(::Type{Variogram}, γ::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares())
  # fit each variogram type
  res = [fit_impl(V, γ, algo) for V in FITTABLE]
  γs, es = first.(res), last.(res)

  # return best candidate
  γs[argmin(es)]
end

function fit_impl(V::Type{<:Variogram}, γ::EmpiricalVariogram, algo::WeightedLeastSquares)
  # values of empirical variogram
  x, y, n = values(γ)

  # custom ball of given radius
  ball(r) = MetricBall(r, distance(γ))

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
  w = f ≠ nothing ? map(f, x) : n / sum(n)

  # objective function
  function J(p)
    g = V(ball(p[1]), sill=p[2] + p[3], nugget=p[3])
    sum(w[i] * (g(x[i]) - y[i])^2 for i in eachindex(x))
  end

  # initial guess
  pₒ = [xmax / 3, 0.95 * ymax, 1e-6]

  # box constraints
  l = [0.0, 0.0, 0.0]
  u = [xmax, ymax, ymax]

  # solve optimization problem
  sol = Optim.optimize(J, l, u, pₒ)
  err = Optim.minimum(sol)
  p = Optim.minimizer(sol)

  # optimal variogram (with units)
  vario = V(ball(p[1]), sill=(p[2] + p[3]) * 𝓊, nugget=p[3] * 𝓊)

  vario, err
end

# convenient method with weighting function as third argument
fit(V, γ::EmpiricalVariogram, weightfun::Function) = fit(V, γ, WeightedLeastSquares(weightfun))
