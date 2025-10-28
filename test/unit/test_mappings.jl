using Plots
using BenchmarkTools
using CurvilinearGrids
using KernelAbstractions

"""
    spow(x, p)

Signed-power
"""
spow(x, p) = sign(x) * abs(x)^p

# function piecewise_mapping(i, N, p, fine_spacing_coef, width, origin, offset)
#   x = -1 + 2 * (i - 1 + offset) / (N - 1)
#   s = width / (1 + fine_spacing_coef)
#   slope = s * (p + fine_spacing_coef)
#   if x < -1
#     y = -width + slope * (x + 1) + origin
#   elseif x > 1
#     y = width + slope * (x - 1) + origin
#   else
#     y = s * (spow(x, p) + fine_spacing_coef * x) + origin
#   end

#   y
# end

"""
    piecewise_mapping(i, N, p, fine_spacing_coef, width, origin, offset)

# Arguments
 - `i`: index
 - `N`: total number of indices
 - `p`: exponent for the inner region
 - `fine_spacing_coef`: slope of the inner region (must be > 0)
 - `width`: approximate width of the inner region
 - `origin`: origin of the inner region
 - `offset`: offset of the wings: >0 biases to the right, <0 biases to the left
"""
function piecewise_smooth(i, N::Int, p, fine_spacing_coef, width, origin, offset, ε)
  # map i ∈ [1, N] → x ∈ [-1, 1]
  x = -1 + 2 * (i - 1 + offset) / (N - 1)

  # normalization
  s = width / (1 + fine_spacing_coef)
  y_c(x) = s * (spow(x, p) + fine_spacing_coef * x)
  dy_c(x) = s * (p * abs(x)^(p - 1) + fine_spacing_coef)

  # values and slopes at boundaries
  y_r, slope_r = y_c(1.0), dy_c(1.0)
  y_l, slope_l = y_c(-1.0), dy_c(-1.0)

  # smooth transition functions
  S_right = 0.5 * (1 + tanh((x - 1) / ε))
  S_left = 0.5 * (1 + tanh((-x - 1) / ε))

  y_right = y_r + slope_r * (x - 1)
  y_left = y_l + slope_l * (x + 1)

  # smooth blending of inner and outer regions
  y =
    (1 - S_right - S_left) * y_c(x) +
    S_right * y_right + #
    S_left * y_left + #
    origin

  return y
end

# x(i) = piecewise_mapping(
#   i,
#   300, # N
#   25, # p
#   0.45, # fine_spacing_coef
#   70.0, # width
#   400.0, # origin
#   2.0, # offset
# )

# x_xi = x.(1:300)

# plot(
#   x_xi,
#   zeros(length(x_xi));
#   marker=:circle,
#   ms=3,
#   size=(800, 200),
#   ylims=(-1, 1),
#   legend=false,
#   yticks=nothing,
# )

function sector_mapping(θmin, θmax, ϕmin, ϕmax, ncells::NTuple{3,Int})
  ni, nj, nk = ncells

  # Uniform spacings
  Δθ = (θmax - θmin) / nj
  Δϕ = (ϕmax - ϕmin) / nk

  # Coordinate functions
  r(i) = piecewise_smooth(
    i,
    ni, # N
    11, # p
    0.45, # fine_spacing_coef
    70.0, # width
    400.0, # origin
    2.0, # offset
    0.05, # ε transition smoothing
  )

  θ(j) = (θmin + (j - 1) * Δθ) / pi
  ϕ(k) = (ϕmin + (k - 1) * Δϕ) / pi

  x(i, j, k) = r(i) * sinpi(θ(j)) * cospi(ϕ(k))
  y(i, j, k) = r(i) * sinpi(θ(j)) * sinpi(ϕ(k))
  z(i, j, k) = r(i) * cospi(θ(j))

  return (x, y, z)
end

celldims = (100, 50, 50)
θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)
ϕmin, ϕmax = -deg2rad(5), deg2rad(5)

x, y, z = sector_mapping(θmin, θmax, ϕmin, ϕmax, celldims);
# backend = AutoForwardDiff()
@profview begin
  mesh = ContinuousCurvilinearGrid3D(x, y, z, celldims, :meg6, CPU())
end
# gcl_identities, max_vals = gcl(mesh.edge_metrics, mesh.iterators.cell.domain, eps());
save_vtk(mesh, "sector_ad_mesh")

gcl_identities, max_vals = gcl(mesh.edge_metrics, mesh.iterators.cell.domain, 5e-13)
