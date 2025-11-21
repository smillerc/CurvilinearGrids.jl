using CurvilinearGrids,
  Test,
  KernelAbstractions,
  DifferentiationInterface,
  BenchmarkTools,
  UnPack,
  CartesianDomains

@testset "Wavy ContinuousCurvilinearGrid3D" begin
  function wavy_mapping()
    function x(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin = p
      xmin +
      Δx0 * ((i - 1) + Ax * sin(pi * n * (j - 1) * Δy0) * sin(pi * n * (k - 1) * Δz0))
    end

    function y(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, ymin = p
      ymin +
      Δy0 * ((j - 1) + Ay * sin(pi * n * (k - 1) * Δz0) * sin(pi * n * (i - 1) * Δx0))
    end

    function z(t, i, j, k, p)
      @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, zmin = p
      zmin +
      Δz0 * ((k - 1) + Az * sin(pi * n * (i - 1) * Δx0) * sin(pi * n * (j - 1) * Δy0))
    end

    return (x, y, z)
  end

  function wavy_params()
    celldims = (41, 41, 41)

    ni, nj, nk = celldims
    Lx = Ly = Lz = 12

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    Ax = 0.2 / Δx0
    Ay = 0.4 / Δy0
    Az = 0.6 / Δz0

    n = 0.5

    params = (; Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin, ymin, zmin)
    return params, celldims
  end

  x, y, z = wavy_mapping()
  params, celldims = wavy_params()
  mesh = ContinuousCurvilinearGrid3D(x, y, z, params, celldims, :meg6, CPU())
  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)
  # @show extrema(I1)
  # @show extrema(I2)
  # @show extrema(I3)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
  @test all(abs.(extrema(I3)) .< 1e-14)
  nothing
end

# begin
#   ξηζ = (
#     (
#       0.0,
#       (mesh.iterators.global_domain.cell.full[I].I .- mesh.iterators.nhalo .+ (1 / 2))...,
#       params,
#     ) for I in mesh.iterators.cell.full
#   )
#   ξηζ

#   mesh.metric_functions_cache.edge.ξ̂xᵢ₊½(first(ξηζ))
# end

@testset "Sphere Sector ContinuousCurvilinearGrid3D" begin
  function get_sector_parameters()
    celldims = (40, 40, 40)
    ni, nj, nk = celldims
    rmin, rmax = 1.0, 4.0
    θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)   # narrow polar band around equator
    ϕmin, ϕmax = -deg2rad(10), deg2rad(10)

    # Uniform spacings
    Δr = (rmax - rmin) / ni
    Δθ = (θmax - θmin) / nj
    Δϕ = (ϕmax - ϕmin) / nk

    params = (; Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin)
    return params, celldims
  end

  function get_sector_mapping()
    function x(t, i, j, k, p)
      @unpack Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin = p
      r(i) = (rmin + (i - 1) * Δr)
      θ(j) = (θmin + (j - 1) * Δθ)
      ϕ(k) = (ϕmin + (k - 1) * Δϕ)
      return r(i) * sin(θ(j)) * cos(ϕ(k))
    end

    function y(t, i, j, k, p)
      @unpack Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin = p
      r(i) = (rmin + (i - 1) * Δr)
      θ(j) = (θmin + (j - 1) * Δθ)
      ϕ(k) = (ϕmin + (k - 1) * Δϕ)
      return r(i) * sin(θ(j)) * sin(ϕ(k))
    end

    function z(t, i, j, k, p)
      @unpack Δr, rmax, rmin, Δθ, θmax, θmin, Δϕ, ϕmax, ϕmin = p
      r(i) = (rmin + (i - 1) * Δr)
      θ(j) = (θmin + (j - 1) * Δθ)
      ϕ(k) = (ϕmin + (k - 1) * Δϕ)
      return r(i) * cos(θ(j))
    end

    return (x, y, z)
  end

  sector_params, celldims = get_sector_parameters()
  x, y, z = get_sector_mapping()

  backend = AutoForwardDiff()
  mesh = ContinuousCurvilinearGrid3D(
    x, y, z, sector_params, celldims, :meg6, CPU(), backend
  )
  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)
  # @show extrema(I1)
  # @show extrema(I2)
  # @show extrema(I3)
  # # CurvilinearGrids.save_vtk(mesh, "spherical_sector_ad")

  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
  @test all(abs.(extrema(I3)) .< 1e-14)
  nothing
end

@testset "Uniform ContinuousCurvilinearGrid3D" begin
  function get_uniform_mapping()
    function x(t, i, j, k, p)
      @unpack xmin, Δx = p
      return xmin + (i - 1) * Δx
    end

    function y(t, i, j, k, p)
      @unpack ymin, Δy = p
      return ymin + (j - 1) * Δy
    end

    function z(t, i, j, k, p)
      @unpack zmin, Δz = p
      return zmin + (k - 1) * Δz
    end

    return (x, y, z, params)
  end

  function get_uniform_params()
    celldims = (40, 80, 120)
    xmin, xmax = (0.0, 2.0)
    ymin, ymax = (1, 3)
    zmin, zmax = (-1, 2)

    ni, nj, nk = celldims

    Δx = (xmax - xmin) / ni
    Δy = (ymax - ymin) / nj
    Δz = (zmax - zmin) / nk

    params = (; xmin, ymin, zmin, Δx, Δy, Δz)
    return params, celldims
  end

  params, celldims = get_uniform_params()
  x, y, z = get_uniform_mapping()

  backend = AutoForwardDiff()
  mesh = ContinuousCurvilinearGrid3D(x, y, z, params, celldims, :meg6, CPU(), backend)
  I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(mesh.edge_metrics, mesh.iterators.cell.domain)
  # @show extrema(I1)
  # @show extrema(I2)
  # @show extrema(I3)
  @test all(abs.(extrema(I1)) .< 1e-14)
  @test all(abs.(extrema(I2)) .< 1e-14)
  @test all(abs.(extrema(I3)) .< 1e-14)

  cell_volume = 0.05 * 0.025 * 0.025

  @test all(mesh.cell_center_metrics.J .≈ cell_volume)

  @test all(mesh.cell_center_metrics.x₁.ξ .≈ 0.05)
  @test all(mesh.cell_center_metrics.x₂.ξ .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₃.ξ .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₁.η .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₂.η .≈ 0.025)
  @test all(mesh.cell_center_metrics.x₃.η .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₁.ζ .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₂.ζ .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₃.ζ .≈ 0.025)

  @test all(mesh.cell_center_metrics.ξ̂.x₁ .≈ 0.000625)
  @test all(mesh.cell_center_metrics.ξ̂.x₂ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ξ̂.x₃ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η̂.x₁ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η̂.x₂ .≈ 0.00125)
  @test all(mesh.cell_center_metrics.η̂.x₃ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ̂.x₁ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ̂.x₂ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ̂.x₃ .≈ 0.00125)

  @test all(mesh.cell_center_metrics.ξ.x₁ .≈ 20.0)
  @test all(mesh.cell_center_metrics.ξ.x₂ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ξ.x₃ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.x₁ .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.x₂ .≈ 40.0)
  @test all(mesh.cell_center_metrics.η.x₃ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₁ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₂ .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₃ .≈ 40.0)
  @test all(mesh.cell_center_metrics.ξ.t .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.t .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.t .≈ 0.0)

  iaxis, jaxis, kaxis = (1, 2, 3)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)
  k₊½_domain = expand(domain, kaxis, -1)

  @test all(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain] .≈ 0.000625)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₁[j₊½_domain] .≈ 0.000625)
  @test all(mesh.edge_metrics.k₊½.ξ̂.x₁[k₊½_domain] .≈ 0.000625)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₂[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₂[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.x₂[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₃[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₃[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.x₃[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.x₁[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.x₁[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.x₁[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.x₂[i₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.j₊½.η̂.x₂[j₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.k₊½.η̂.x₂[k₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.i₊½.η̂.x₃[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.x₃[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.x₃[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.x₁[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.x₁[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.x₁[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.x₂[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.x₂[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.x₂[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.x₃[i₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.j₊½.ζ̂.x₃[j₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.k₊½.ζ̂.x₃[k₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.t[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.t[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.t[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.t[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.t[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.t[k₊½_domain] .≈ 0.0)
end

# @testset "ContinuousCurvilinearGrid3D vs CurvilinearGrid3D" begin
#   nhalo = 5
#   nr, nθ, nϕ = 50, 51, 51
#   celldims = (nr, nθ, nϕ)

#   rmin, rmax = (0.0400, 0.05500)
#   θ0 = 90.0
#   ϕ0 = 0.0

#   Δθ = 20.0
#   Δϕ = 20.0

#   θmin, θmax = ((θ0 + Δθ, θ0 - Δθ) .|> deg2rad)
#   ϕmin, ϕmax = ((ϕ0 + Δϕ, ϕ0 - Δϕ) .|> deg2rad)

#   # rmin, rmax = 1.0, 4.0
#   # θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)   # narrow polar band around equator
#   # ϕmin, ϕmax = -deg2rad(10), deg2rad(10)

#   (x, y, z) = spherical_sector_mapping(rmin, rmax, θmin, θmax, ϕmin, ϕmax, celldims)
#   # (x, y, z) = wavy_mapping(celldims)

#   # diff_backend = AutoForwardDiff()
#   @info "ContinuousCurvilinearGrid3D"
#   cm = ContinuousCurvilinearGrid3D(x, y, z, celldims, :meg6, CPU())

#   I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(cm.edge_metrics, cm.iterators.cell.domain)
#   # @show extrema(I1)
#   # @show extrema(I2)
#   # @show extrema(I3)
#   @test all(abs.(extrema(I1)) .< 1e-19)
#   @test all(abs.(extrema(I2)) .< 1e-19)
#   @test all(abs.(extrema(I3)) .< 5e-14)

#   # CurvilinearGrids.save_vtk(cm, "spherical_sector_ad")

#   @info "CurvilinearGrid3D"
#   xdom = cm.node_coordinates.x[cm.iterators.node.full]
#   ydom = cm.node_coordinates.y[cm.iterators.node.full]
#   zdom = cm.node_coordinates.z[cm.iterators.node.full]
#   dm = CurvilinearGrid3D(xdom, ydom, zdom, :meg6; halo_coords_included=true)
#   I1, I2, I3 = CurvilinearGrids.GridTypes.gcl(dm.edge_metrics, dm.iterators.cell.domain)
#   # @show extrema(I1)
#   # @show extrema(I2)
#   # @show extrema(I3)
#   @test all(abs.(extrema(I1)) .< 5e-14)
#   @test all(abs.(extrema(I2)) .< 5e-14)
#   @test all(abs.(extrema(I3)) .< 5e-14)

#   # CurvilinearGrids.save_vtk(cm, "spherical_sector_meg")

#   dom = cm.iterators.cell.domain

#   # for dim in (:ξ, :η, :ζ, :ξ̂, :η̂, :ζ̂, :x₁, :x₂, :x₃)
#   #   for ((dm_name, dm_component), (cm_name, cm_component)) in zip(
#   #     StructArrays.components(dm.cell_center_metrics[dim]) |> pairs,
#   #     StructArrays.components(cm.cell_center_metrics[dim]) |> pairs,
#   #   )
#   #     # @test all(isapprox.(dm_component[dom], cm_component[dom], rtol=1e-5))

#   #     passes = all(isapprox.(dm_component[dom], cm_component[dom], rtol=1e-3))
#   #     @info "Dim: $dim, $dm_name, passes? $passes"
#   #     if !passes
#   #       @show extrema(dm_component[dom])
#   #       @show extrema(cm_component[dom])
#   #       println()
#   #     end
#   #   end
#   #   # println()
#   # end

#   nothing
# end
