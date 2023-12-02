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
    fit(V, g, [algo])
    fit(V, g, [weightfun])

Fit theoretical variogram type `V` to empirical variogram `g`
using algorithm `algo`. Default algorithm is `WeightedLeastSquares`.

Alternatively pass the weighting function `weightfun` directly
to the fitting procedure.
"""
fit(V::Type{<:Variogram}, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares()) = fit_impl(V, g, algo) |> first

"""
    fit(Variogram, g, [algo])
    fit(Variogram, g, [weightfun])

Fit all subtypes of `Variogram` to empirical variogram `g` and
return the one with minimum error as defined by the algorithm `algo`.

Alternatively pass the weighting function `weightfun` directly
to the fitting procedure.
"""
function fit(::Type{Variogram}, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares())
  # fit each variogram type
  res = [fit_impl(V, g, algo) for V in FITTABLE]
  γs, ϵs = first.(res), last.(res)

  # return best candidate
  γs[argmin(ϵs)]
end

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

# convenient method with weighting function as third argument
fit(V, g::EmpiricalVariogram, weightfun::Function) = fit(V, g, WeightedLeastSquares(weightfun))
