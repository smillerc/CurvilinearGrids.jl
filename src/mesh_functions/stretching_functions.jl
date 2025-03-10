"""
    one_sided_stretch((x0, x1), N::Int, α)

Apply a one-sided tanh stretching function. The `(x0, x1)` points
represent the desired end-points, `N` is the number of points, and
`α` is the stretching factor (recommended to be no more than ~10).
The default option is to cluster near `x0`, but setting `cluster_on_x0=false`
will cluster near `x1`
"""
function one_sided_stretch((x0, x1), N, α; cluster_on_x0=true)
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
    double_sided_stretch((x0, x1), N::Int, α)

Apply a doubled-sided stretching function, where points are clustered
on either endpoint. The stretching factor `α` must be between 0 and 1
"""
function double_sided_stretch((x0, x1), N, α)
  @assert 0 <= α <= 1 "The stretching factor α must be between 0 and 1"
  len = abs(x1 - x0)

  s(ξ) = (1 + tanh(α * (ξ / N - (1 / 2))) / (tanh(α / 2))) / 2
  p = s.(0:(N - 1))
  return p .* len .+ x0
end