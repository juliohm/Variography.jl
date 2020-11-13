# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(v::EmpiricalVarioplane; showrange=true, model=GaussianVariogram)
  # retrieve info
  θs = v.θs
  γs = v.γs

  # abscissa is the same for all variograms
  rs, _, __ = values(γs[1])

  # variogram values for all variograms
  zs = map(γs) do γ
    _, zs, __ = values(γ)

    # handle NaN values (i.e. empty bins)
    isnan(zs[1]) && (zs[1] = 0)
    for i in 2:length(zs)
      isnan(zs[i]) && (zs[i] = zs[i-1])
    end

    zs
  end
  Z = reduce(hcat, zs)

  # exploit symmetry
  θs = range(0, 2π, length=2*length(θs))
  Z  = [Z Z]

  # plot settings
  projection --> :polar
  ylims --> (0, Inf)

  # plot varioplane
  @series begin
    seriestype --> :heatmap
    θs, rs, Z
  end

  # plot ranges
  if showrange
    ls = [range(fit(model, γ)) for γ in γs]
    @series begin
      seriestype --> :path
      linecolor --> :black
      label --> "range"
      θs, [ls; ls]
    end
  end
end
