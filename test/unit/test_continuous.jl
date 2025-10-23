using CurvilinearGrids
using Test
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

@testset "ContinuousCurvilinearGrid3D vs CurvilinearGrid3D" begin
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
  cm = ContinuousCurvilinearGrid3D(x, y, z, celldims, :meg6, CPU(), AutoForwardDiff())

  @info "CurvilinearGrid3D"
  xdom = cm.node_coordinates.x[cm.iterators.node.domain]
  ydom = cm.node_coordinates.y[cm.iterators.node.domain]
  zdom = cm.node_coordinates.z[cm.iterators.node.domain]
  dm = CurvilinearGrid3D(xdom, ydom, zdom, :meg6)

  # save_vtk(dm, "sector_meg")
  # save_vtk(cm, "sector_ad")

  gcl(cm.edge_metrics, cm.iterators.cell.domain)
  gcl(dm.edge_metrics, dm.iterators.cell.domain)

  for dim in (:ξ, :η, :ζ, :ξ̂, :η̂, :ζ̂, :x₁, :x₂, :x₃)
    for ((dm_name, dm_component), (cm_name, cm_component)) in zip(
      StructArrays.components(dm.cell_center_metrics[dim]) |> pairs,
      StructArrays.components(cm.cell_center_metrics[dim]) |> pairs,
    )
      @test all(dm_component[dom] .≈ cm_component[dom])
      @info "Dim: $dim, $dm_name, passes? $passes"
      if !passes
        @show extrema(dm_component[dom])
        @show extrema(cm_component[dom])
        println()
      end
    end
    println()
  end

  nothing
end
