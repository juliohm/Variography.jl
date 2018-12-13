# ------------------------------------------------------------------
# Copyright (c) 2017, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

@userplot VarPlane

@recipe function f(vp::VarPlane; normal=nothing,
                   nlags=20, maxlag=nothing,
                   nangs=50, atol=10., btol=0.95,
                   showrange=true, varmodel=GaussianVariogram)
  # get inputs
  spatialdata = vp.args[1]
  var₁ = vp.args[2]
  var₂ = length(vp.args) == 3 ? vp.args[3] : var₁

  # polar coordinates and variogram values
  rs = Vector{Float64}(undef, nlags)
  θs = range(0, stop=π, length=nangs)
  zs = Matrix{Float64}(undef, nlags, nangs)

  # ranges for each angle
  ls = Vector{Float64}(undef, nangs)

  # loop over half plane
  for (j, θ) in Iterators.enumerate(θs)
    γ = DirectionalVariogram(spatialdata, (cos(θ),sin(θ)), var₁, var₂;
                             nlags=nlags, maxlag=maxlag, atol=atol, btol=btol)
    rs, ys, _ = values(γ)

    # clean NaN values (i.e. empty bins)
    isnan(ys[1]) && (ys[1] = 0)
    for i in 2:length(ys)
      isnan(ys[i]) && (ys[i] = ys[i-1])
    end

    # save variogram for this angle
    zs[:,j] = ys

    if showrange
      γtheo = fit(varmodel, γ)
      ls[j] = range(γtheo)
    end
  end

  # exploit symmetry
  θs = range(0, stop=2π, length=2nangs)
  zs = hcat(zs, zs)

  projection --> :polar
  ylim --> (0, Inf)

  # plot variogram plane
  @series begin
    seriestype --> :heatmap
    color --> :curl

    θs, rs, zs
  end

  # plot ranges
  if showrange
    @series begin
      seriestype --> :path
      linecolor --> :black
      label --> "range"

      θs, ls
    end
  end
end
