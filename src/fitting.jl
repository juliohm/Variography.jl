# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    FitAlgo

An algorithm for fitting theoretical variograms.
"""
abstract type FitAlgo end

"""
    WeightedLeastSquares(weightfun)

Fit theoretical variogram using weighted least squares.
The default weighting function `weightfun` is `x -> 1`
(i.e. no weighting)
"""
struct WeightedLeastSquares <: FitAlgo
  weightfun::Function
end

WeightedLeastSquares() = WeightedLeastSquares(x -> 1.)

"""
    fit(V, γ, [algo])

Fit theoretical variogram type `V` to empirical variogram `γ`
using algorithm `algo`. Default algorithm is `WeightedLeastSquares`.
"""
function fit(::Type{V}, γ::EmpiricalVariogram,
             algo::FitAlgo=WeightedLeastSquares()) where {V<:Variogram}
  # dispatch appropriate implementation
  vario, err = fit_impl(V, γ, algo)

  vario
end

"""
    fit(Variogram, γ, [algo])

Fit all stationary variogram types to empirical variogram `γ`,
which are subtypes of `Variogram`, and return the one with
minimum error as defined by the algorithm `algo`.
"""
function fit(::Type{Variogram}, γ::EmpiricalVariogram,
             algo::FitAlgo=WeightedLeastSquares())
  # list of variogram types to try
  Vs = (GaussianVariogram, SphericalVariogram, ExponentialVariogram,
        MaternVariogram, CubicVariogram, PentasphericalVariogram,
        SineHoleVariogram)

  # fit each variogram type
  res    = [fit_impl(V, γ, algo) for V in Vs]
  varios = first.(res)
  errs   = last.(res)

  # return best candidate
  varios[argmin(errs)]
end

function fit_impl(::Type{V}, γ::EmpiricalVariogram,
                  algo::WeightedLeastSquares) where {V<:Variogram}
  # retrieve empirical points
  x, y, n = values(γ)

  # discard invalid bins
  x = x[n .> 0]
  y = y[n .> 0]

  # evaluate weights
  w = map(algo.weightfun, x)

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
  sol = optimize(J, l, u, pₒ)
  err = Optim.minimum(sol)
  p   = Optim.minimizer(sol)

  # best variogram
  vario = V(range=p[1], sill=p[2], nugget=p[3])

  vario, err
end
