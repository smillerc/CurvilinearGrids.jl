
function uniform_mapping(xmin, xmax, ymin, ymax, zmin, zmax, ncells::NTuple{3,Int})
  ni, nj, nk = ncells

  Δx = (xmax - xmin) / ni
  Δy = (ymax - ymin) / nj
  Δz = (zmax - zmin) / nk

  x(i, j, k) = xmin + (i - 1) * Δx
  y(i, j, k) = ymin + (j - 1) * Δy
  z(i, j, k) = zmin + (k - 1) * Δz

  return (x, y, z)
end

function spherical_sector_mapping(rmin, rmax, θmin, θmax, ϕmin, ϕmax, ncells::NTuple{3,Int})
  ni, nj, nk = ncells

  # Uniform spacings
  Δr = (rmax - rmin) / ni
  Δθ = (θmax - θmin) / nj
  Δϕ = (ϕmax - ϕmin) / nk

  # Coordinate functions
  r(i) = rmin + (i - 1) * Δr
  θ(j) = (θmin + (j - 1) * Δθ) / pi
  ϕ(k) = (ϕmin + (k - 1) * Δϕ) / pi

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

  Ax = 0.0 / Δx0
  Ay = 0.0 / Δy0
  Az = 0.0 / Δz0

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

@testset "Wavy ContinuousCurvilinearGrid3D" begin
  celldims = (41, 41, 41)

  (x, y, z) = wavy_mapping(celldims)

  backend = AutoForwardDiff()
  mesh = ContinuousCurvilinearGrid3D(x, y, z, celldims, :meg6, CPU())
  gcl_identities, max_vals = gcl(mesh.edge_metrics, mesh.iterators.cell.domain, eps())
  # @show gcl_identities, max_vals
  @test all(gcl_identities)
end

@testset "Sphere Sector ContinuousCurvilinearGrid3D" begin
  nhalo = 5
  nr, nθ, nϕ = 41, 41, 41
  celldims = (nr, nθ, nϕ)
  rmin, rmax = 1.0, 4.0
  θmin, θmax = π / 2 - deg2rad(5), π / 2 + deg2rad(5)   # narrow polar band around equator
  ϕmin, ϕmax = -deg2rad(10), deg2rad(10)

  (x, y, z) = spherical_sector_mapping(rmin, rmax, θmin, θmax, ϕmin, ϕmax, celldims)

  backend = AutoForwardDiff()
  mesh = ContinuousCurvilinearGrid3D(x, y, z, celldims, :meg6, CPU())
  gcl_identities, max_vals = gcl(mesh.edge_metrics, mesh.iterators.cell.domain, eps())
  # @show gcl_identities, max_vals
  @test all(gcl_identities)
end

@testset "Uniform ContinuousCurvilinearGrid3D" begin
  x0, x1 = (0.0, 2.0)
  y0, y1 = (1, 3)
  z0, z1 = (-1, 2)
  celldims = (40, 80, 120)
  (x, y, z) = uniform_mapping(x0, x1, y0, y1, z0, z1, celldims)

  backend = AutoForwardDiff()
  mesh = ContinuousCurvilinearGrid3D(x, y, z, celldims, :meg6, CPU())
  gcl_identities, max_vals = gcl(mesh.edge_metrics, mesh.iterators.cell.domain, eps())
  # @show gcl_identities, max_vals
  @test all(gcl_identities)

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
  cm = ContinuousCurvilinearGrid3D(x, y, z, celldims, :meg6, CPU())

  @info "CurvilinearGrid3D"
  xdom = cm.node_coordinates.x[cm.iterators.node.full]
  ydom = cm.node_coordinates.y[cm.iterators.node.full]
  zdom = cm.node_coordinates.z[cm.iterators.node.full]
  dm = CurvilinearGrid3D(xdom, ydom, zdom, :meg6; halo_coords_included=true)

  # save_vtk(dm, "sector_meg")
  # save_vtk(cm, "sector_ad")

  cm_gcl_identities, cm_max_vals = gcl(cm.edge_metrics, cm.iterators.cell.domain, eps())
  dm_gcl_identities, dm_max_vals = gcl(dm.edge_metrics, dm.iterators.cell.domain, eps())
  @test all(cm_gcl_identities)
  @test all(dm_gcl_identities)

  dom = cm.iterators.cell.domain

  for dim in (:ξ, :η, :ζ, :ξ̂, :η̂, :ζ̂, :x₁, :x₂, :x₃)
    for ((dm_name, dm_component), (cm_name, cm_component)) in zip(
      StructArrays.components(dm.cell_center_metrics[dim]) |> pairs,
      StructArrays.components(cm.cell_center_metrics[dim]) |> pairs,
    )
      @test all(dm_component[dom] .≈ cm_component[dom])
      # passes = all(dm_component[dom] .≈ cm_component[dom])
      # @info "Dim: $dim, $dm_name, passes? $passes"
      # if !passes
      #   @show extrema(dm_component[dom])
      #   @show extrema(cm_component[dom])
      #   println()
      # end
    end
    # println()
  end

  nothing
end
