using Roots

"""
    one_sided_stretch((x0, x1), N::Int, α)

Apply a one-sided tanh stretching function. The `(x0, x1)` points
represent the desired end-points, `N` is the number of points, and
`α` is the stretching factor (recommended to be no more than ~10).
The default option is to cluster near `x0`, but setting `cluster_on_x0=false`
will cluster near `x1`
"""
function one_sided_stretch((x0, x1), N::Int, α; cluster_on_x0=true)
  len = abs(x1 - x0)
  s(ξ) = 1 + tanh(α * ((ξ / (N - 1)) - 1)) / tanh(α)

  if cluster_on_x0
    p = s.(0:(N - 1))
  else
    p = reverse(s.(0:(N - 1)))
  end
  return p .* len .+ x0
end

"""
    one_sided_with_initial_spacing((x0, x1), N::Int, dx0; cluster_on_x0=true)

Apply a one-sided tanh stretching function but specify the spacing at the first layer.
The `(x0, x1)` points represent the desired end-points, `N` is the number of points.
The default option is to cluster near `x0`, but setting `cluster_on_x0=false`
will cluster near `x1`.
"""
function one_sided_with_initial_spacing((x0, x1), N::Int, dx0; cluster_on_x0=true)
  if dx0 >= abs(x1 - x0) / N
    @warn "Desired spacing dx0 = $dx0 is >= uniform spacing (abs(x1 - x0) / N). Just use a uniform range instead of the one-sided stretch function."
    return nothing
  end

  function g(α)
    x = one_sided_stretch((x0, x1), N, α; cluster_on_x0=cluster_on_x0)
    dx = x[2] - x[1]
    objective = dx0 - dx
    return objective
  end

  α_needed = find_zero(g, 1.0)

  return one_sided_stretch((x0, x1), N, α_needed; cluster_on_x0=cluster_on_x0)
end

"""
    double_sided_stretch((x0, x1), N::Int, α)

Apply a doubled-sided stretching functioN::Int, where points are clustered
on either endpoint. The stretching factor `α` must be between 0 and 1
"""
function double_sided_stretch((x0, x1), N::Int, α)
  @assert 0 <= α <= 1 "The stretching factor α must be between 0 and 1"
  len = abs(x1 - x0)

  s(ξ) = (1 + tanh(α * (ξ / N - (1 / 2))) / (tanh(α / 2))) / 2
  p = s.(0:(N - 1))
  return p .* len .+ x0
end
