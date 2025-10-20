using CurvilinearGrids
using DifferentiationInterface
using ForwardDiff: ForwardDiff
using BenchmarkTools
using StaticArrays
using KernelAbstractions

function spherical_sector_mapping(
  rmin, rmax, θmin, θmax, ϕmin, ϕmax, ncells::NTuple{3,Int}, nhalo
)
  ni, nj, nk = ncells

  # Uniform spacings
  Δr = (rmax - rmin) / ni
  Δθ = (θmax - θmin) / nj
  Δϕ = (ϕmax - ϕmin) / nk

  # Coordinate functions
  r(i) = rmin + (i - 1) * Δr
  θ(j) = θmin + (j - 1) * Δθ
  ϕ(k) = ϕmin + (k - 1) * Δϕ

  x(i, j, k) = r(i) * sin(θ(j)) * cos(ϕ(k))
  y(i, j, k) = r(i) * sin(θ(j)) * sin(ϕ(k))
  z(i, j, k) = r(i) * cos(θ(j))

  return (x, y, z)
end

function wavy_mapping(ncells::NTuple{3,Int})
  ni, nj, nk = ncells
  Lx = Ly = Lz = 12

  xmin = -Lx / 2
  ymin = -Ly / 2
  zmin = -Lz / 2

  Δx0 = Lx / ni
  Δy0 = Ly / nj
  Δz0 = Lz / nk

  Ax = 0.2 / Δx0
  Ay = 0.2 / Δy0
  Az = 0.2 / Δz0

  n = 0.5

  function x(i, j, k)
    xmin + Δx0 * ((i - 1) + Ax * sin(pi * n * (j - 1) * Δy0) * sin(pi * n * (k - 1) * Δz0))
  end
  function y(i, j, k)
    ymin + Δy0 * ((j - 1) + Ay * sin(pi * n * (k - 1) * Δz0) * sin(pi * n * (i - 1) * Δx0))
  end
  function z(i, j, k)
    zmin + Δz0 * ((k - 1) + Az * sin(pi * n * (i - 1) * Δx0) * sin(pi * n * (j - 1) * Δy0))
  end

  return (x, y, z)
end

function gcl(
  em, # edge metrics
  domain,
  ϵ=5e-13,
)
  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  I₁ = -Inf
  I₂ = -Inf
  I₃ = -Inf
  for idx in domain
    i, j, k = idx.I
    _I₁ = (
      (em.i₊½.ξ̂.x₁[i, j, k] - em.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (em.j₊½.η̂.x₁[i, j, k] - em.j₊½.η̂.x₁[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₁[i, j, k] - em.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    _I₂ = (
      (em.i₊½.ξ̂.x₂[i, j, k] - em.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (em.j₊½.η̂.x₂[i, j, k] - em.j₊½.η̂.x₂[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₂[i, j, k] - em.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    _I₃ = (
      (em.i₊½.ξ̂.x₃[i, j, k] - em.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (em.j₊½.η̂.x₃[i, j, k] - em.j₊½.η̂.x₃[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₃[i, j, k] - em.k₊½.ζ̂.x₃[i, j, k - 1])
    )

    I₁_passes = abs(_I₁) < ϵ
    I₂_passes = abs(_I₂) < ϵ
    I₃_passes = abs(_I₃) < ϵ

    I₁ = max(I₁, abs(_I₁))
    I₂ = max(I₂, abs(_I₂))
    I₃ = max(I₃, abs(_I₃))
    # @show idx, I₁, I₂, I₃
    # break
    if !(I₁_passes && I₂_passes && I₃_passes)
      break
    end
  end
  @show I₁ I₂ I₃

  return (I₁_passes, I₂_passes, I₃_passes)
end

begin
  nhalo = 5
  nr, nθ, nϕ = 20, 20, 20
  celldims = (nr, nθ, nϕ)
  rmin, rmax = 1.0, 4.0
  θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)   # narrow polar band around equator
  ϕmin, ϕmax = -deg2rad(10), deg2rad(10)

  (x, y, z) = spherical_sector_mapping(rmin, rmax, θmin, θmax, ϕmin, ϕmax, celldims, nhalo)
  #   (x, y, z) = wavy_mapping(celldims)

  backend = AutoForwardDiff()
  @info "ContinuousCurvilinearGrid3D"
  continuous_mesh = ContinuousCurvilinearGrid3D(
    x, y, z, celldims, nhalo, CPU(), AutoForwardDiff()
  )

  @info "CurvilinearGrid3D"
  discrete_mesh = CurvilinearGrid3D(
    continuous_mesh.node_coordinates..., :meg6; halo_coords_included=true
  )

  #   prep = prepare_derivative(x, backend, 0, Constant(1), Constant(1))
  #   y_grad_prep = prepare_gradient(y, backend, [1, 2, 3])

  #   xξ_1(i, j, k) = ForwardDiff.derivative(ξ -> x(ξ, j, k), i)
  #   xξ_2(i, j, k) = derivative(x, backend, i, Constant(j), Constant(k))
  #   xξ_3(i, j, k) = derivative(x, prep, backend, i, Constant(j), Constant(k))

  #   function xξ_4(i, j, k)
  #     ijk = @SVector [i, j, k]
  #     grad = gradient(x, backend, ijk)
  #     return grad[1]
  #   end

  #   y_ζ(k, i, j) = y(i, j, k)
  #   y_deriv_prep = prepare_derivative(y_ζ, backend, 0, Constant(1), Constant(1))
  #   yζ_1(i, j, k) = ForwardDiff.derivative(ζ -> y(i, j, ζ), k)
  #   yζ_2(i, j, k) = derivative(ζ -> y(i, j, ζ), backend, k)
  #   #   yζ_3(i, j, k) = derivative(y, backend, k, Constant(i), Constant(j))

  #   yζ_3(i, j, k) = derivative(y_ζ, y_deriv_prep, backend, k, Constant(i), Constant(j))
  #   yζ_4(i, j, k) = derivative(y_ζ, AutoForwardDiff(), k, Constant(i), Constant(j))

  #   #   function yζ_3(i, j, k)
  #   #     ijk = [i, j, k]
  #   #     grad = gradient(y, backend, ijk)
  #   #     return grad[3]
  #   #   end

  #   i1, j1, k1 = (10, 15, 18)
  #   #   @show xξ_1(i1, j1, k1)
  #   #   @show xξ_2(i1, j1, k1)
  #   #   @show xξ_3(i1, j1, k1)

  #   @show yζ_1(i1, j1, k1)
  #   @show yζ_2(i1, j1, k1)
  #   @show yζ_3(i1, j1, k1)
  #   @show yζ_4(i1, j1, k1)
  #   @show yζ_5(i1, j1, k1)
  nothing
end

# @benchmark xξ_1($i1, $j1, $k1)
# @benchmark xξ_2($i1, $j1, $k1)
# @benchmark xξ_3($i1, $j1, $k1)
# @benchmark xξ_4($i1, $j1, $k1)

# @benchmark yζ_1($i1, $j1, $k1)
# @benchmark yζ_2($i1, $j1, $k1)
# @benchmark yζ_3($i1, $j1, $k1)
# @benchmark yζ_4($i1, $j1, $k1)

# @code_warntype yζ_4(i1, j1, k1)

gcl(continuous_mesh.metrics.edge, continuous_mesh.iterators.cell.domain)
gcl(discrete_mesh.edge_metrics, discrete_mesh.iterators.cell.domain)

dom = continuous_mesh.iterators.cell.domain
extrema(
  continuous_mesh.metrics.cell.ξ̂.x₂[dom] .- discrete_mesh.cell_center_metrics.ξ̂.x₂[dom]
)
extrema(
  continuous_mesh.metrics.cell.ξ̂.x₃[dom] .- discrete_mesh.cell_center_metrics.ξ̂.x₃[dom]
)
extrema(
  continuous_mesh.metrics.cell.ζ̂.x₃[dom] .- discrete_mesh.cell_center_metrics.ζ̂.x₃[dom]
)
