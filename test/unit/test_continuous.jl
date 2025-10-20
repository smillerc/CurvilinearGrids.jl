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
  θ(j) = (θmin + (j - 1) * Δθ) / pi
  ϕ(k) = (ϕmin + (k - 1) * Δϕ) / pi

  x(i, j, k) = r(i) * sinpi(θ(j)) * cospi(ϕ(k))
  y(i, j, k) = r(i) * sinpi(θ(j)) * sinpi(ϕ(k))
  z(i, j, k) = r(i) * cospi(θ(j))

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
  nr, nθ, nϕ = 21, 21, 21
  celldims = (nr, nθ, nϕ)
  rmin, rmax = 1.0, 4.0
  θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)   # narrow polar band around equator
  ϕmin, ϕmax = -deg2rad(10), deg2rad(10)

  # (x, y, z) = spherical_sector_mapping(rmin, rmax, θmin, θmax, ϕmin, ϕmax, celldims, nhalo)
  (x, y, z) = wavy_mapping(celldims)

  backend = AutoForwardDiff()
  @info "ContinuousCurvilinearGrid3D"
  cm = ContinuousCurvilinearGrid3D(x, y, z, celldims, nhalo, CPU(), AutoForwardDiff())

  @info "CurvilinearGrid3D"
  dm = CurvilinearGrid3D(cm.node_coordinates..., :meg6; halo_coords_included=true)

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

  save_vtk(dm, "sector_meg")
  save_vtk(cm, "sector_ad")
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

gcl(cm.edge_metrics, cm.iterators.cell.domain)
gcl(dm.edge_metrics, dm.iterators.cell.domain)

dom = cm.iterators.cell.domain
extrema(cm.cell_center_metrics.ξ̂.x₂[dom] .- dm.cell_center_metrics.ξ̂.x₂[dom])
extrema(cm.cell_center_metrics.ξ̂.x₃[dom] .- dm.cell_center_metrics.ξ̂.x₃[dom])
extrema(cm.cell_center_metrics.ζ̂.x₃[dom] .- dm.cell_center_metrics.ζ̂.x₃[dom])

cm.metrics.cell.ξ̂.x₂[16, 16, 16]

begin
  i, j, k = (8, 16, 16)
  cm.metrics.edge.j₊½.ξ̂.x₂[i, j, k]

  (
    cm.centroid_coordinates.x[i, j, k],
    cm.centroid_coordinates.y[i, j, k],
    cm.centroid_coordinates.z[i, j, k],
  )
end

cm.metrics.edge.j₊½.ζ̂.x₁[i, j, k] - cm.metrics.edge.j₊½.ζ̂.x₁[i, j - 1, k]
cm.metrics.edge.j₊½.ζ̂.x₂[i, j, k] - cm.metrics.edge.j₊½.ζ̂.x₂[i, j - 1, k]
cm.metrics.edge.j₊½.ζ̂.x₃[i, j, k] - cm.metrics.edge.j₊½.ζ̂.x₃[i, j - 1, k]
cm.metrics.edge.k₊½.ζ̂.x₁[i, j, k] - cm.metrics.edge.k₊½.ζ̂.x₁[i, j, k - 1]
cm.metrics.edge.k₊½.ζ̂.x₂[i, j, k] - cm.metrics.edge.k₊½.ζ̂.x₂[i, j, k - 1]
cm.metrics.edge.k₊½.ζ̂.x₃[i, j, k] - cm.metrics.edge.k₊½.ζ̂.x₃[i, j, k - 1]

dm.edge_metrics.j₊½.ζ̂.x₁[i, j, k] - dm.edge_metrics.j₊½.ζ̂.x₁[i, j - 1, k]
dm.edge_metrics.j₊½.ζ̂.x₂[i, j, k] - dm.edge_metrics.j₊½.ζ̂.x₂[i, j - 1, k]
dm.edge_metrics.j₊½.ζ̂.x₃[i, j, k] - dm.edge_metrics.j₊½.ζ̂.x₃[i, j - 1, k]

dm.edge_metrics.k₊½.ζ̂.x₁[i, j, k] - dm.edge_metrics.k₊½.ζ̂.x₁[i, j, k - 1]
dm.edge_metrics.k₊½.ζ̂.x₂[i, j, k] - dm.edge_metrics.k₊½.ζ̂.x₂[i, j, k - 1]
dm.edge_metrics.k₊½.ζ̂.x₃[i, j, k] - dm.edge_metrics.k₊½.ζ̂.x₃[i, j, k - 1]

# using SpecialFunctions  # for erf

# begin
#   function gaussian(x, x0, fwhm, p)
#     σ = fwhm / (2sqrt(2log(2)))
#     σ² = σ * σ
#     xc = (x - x0)^2
#     return exp(-(((xc / (2σ²)))^p))
#   end

#   # function stretched_r(i, N, rmin, rmax, r0; A=0.2, w=0.1, n=4)
#   #   ξ = (i - 1) / (N - 1)
#   #   ξ₀ = (r0 - rmin) / (rmax - rmin)

#   #   # smooth, continuous warping function
#   #   # warp(ξ) = ξ - A * 0.5 * (erf(((ξ - ξ₀) / w)) + 1)
#   #   # warp(ξ) = ξ - A * (1 / (1 + exp(((ξ - ξ₀) / w)^(2n))))
#   #   warp(ξ) = gaussian(ξ, ξ₀, w, n)

#   #   fξ = warp(ξ)
#   #   return rmin + (rmax - rmin) * fξ
#   # end

#   function stretched_r(i, N, rmin, rmax, r0; A=0.3, fwhm=0.1, p=4)
#     ξ = (i - 1) / (N - 1)
#     ξ0 = (r0 - rmin) / (rmax - rmin)

#     function f(ξ)
#       ξ - A * gaussian(ξ, ξ0, fwhm, p)
#     end
#     f0, f1 = f(0), f(1)
#     Fξ = (f(ξ) - f0) / (f1 - f0)

#     return rmin + (rmax - rmin) * Fξ
#   end

#   # function stretched_r(i, N, rmin, rmax, r0; A=0.3, w=0.1, n=4)
#   #   ξ = (i - 1) / (N - 1)
#   #   ξ₀ = (r0 - rmin) / (rmax - rmin)

#   #   f(ξ) = ξ - A * exp(-((ξ - ξ₀) / w)^(2n))
#   #   f0, f1 = f(0), f(1)
#   #   Fξ = (f(ξ) - f0) / (f1 - f0)

#   #   return rmin + (rmax - rmin) * Fξ
#   # end

#   # function stretched_r_gaussian(
#   #   i::Int, N::Int, rmin::Real, rmax::Real, r0::Real; A::Real=0.4, w::Real=0.06
#   # )
#   #   ξ = (i - 1) / (N - 1)
#   #   ξ0 = (r0 - rmin) / (rmax - rmin)

#   #   # analytic antiderivative terms
#   #   constC = (w * sqrt(pi) / 2)  # convenience
#   #   erf_from0 = erf(-ξ0 / w)     # erf((0-ξ0)/w)

#   #   fξ = ξ - A * constC * (erf((ξ - ξ0) / w) - erf_from0)
#   #   f0 = 0.0 - A * constC * (erf_from0 - erf_from0)   # simplifies to 0.0
#   #   # but compute f1 exactly
#   #   f1 = 1 - A * constC * (erf((1 - ξ0) / w) - erf_from0)

#   #   Fξ = (fξ - f0) / (f1 - f0)
#   #   return rmin + (rmax - rmin) * Fξ
#   # end

#   # function stretched_r_sech2(i, N, rmin, rmax, r0; A::Real=0.5, w::Real=0.05)
#   #   # @assert 1 <= i <= N
#   #   # @assert 0 <= A < 1 "A must be in [0,1) for monotonicity"
#   #   # @assert w > 0
#   #   ξ = (i - 1) / (N - 1)
#   #   ξ0 = (r0 - rmin) / (rmax - rmin)

#   #   t0 = tanh(-ξ0 / w)
#   #   gξ = ξ - A * w * (tanh((ξ - ξ0) / w) - t0)   # analytic antiderivative with C chosen so g(0)=0
#   #   g1 = 1 - A * w * (tanh((1 - ξ0) / w) - t0)
#   #   Gξ = gξ / g1

#   #   return rmin + (rmax - rmin) * Gξ
#   # end

#   rmin, rmax = 0.0, 10.0
#   r0 = 3.5
#   N = 100

#   r = [stretched_r(i, N, rmin, rmax, r0; A=0.2, fwhm=0.5, p=4) for i in 1:N]
#   # r = [stretched_r_sech2(i, N, rmin, rmax, r0; A=0.2, w=0.095) for i in 1:N]
#   # r = [stretched_r_gaussian(i, N, rmin, rmax, r0; A=0.2, w=0.1) for i in 1:N]
#   plot(diff(r); marker=:square)
#   # plot(r, zeros(size(r)); marker=:square)
# end