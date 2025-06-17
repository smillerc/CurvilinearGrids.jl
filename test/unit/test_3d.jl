
# @testset "3D Rectangular Mesh Metrics, Conserved Metrics, GCL"
# using CurvilinearGrids
# begin

#   x0, x1 = (0.0, 2.0)
#   y0, y1 = (1.0, 3.0)
#   z0, z1 = (-1.0, 2.0)
#   ni, nj, nk = (40, 80, 120)

#   mesh = rectilinear_grid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :meg6)
#   @code_warntype rectilinear_grid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :meg6)

#   nothing
# end

@testset "3D Rectangular Mesh Metrics, Conserved Metrics, GCL" begin
  x0, x1 = (0.0, 2.0)
  y0, y1 = (1, 3)
  z0, z1 = (-1, 2)
  ni, nj, nk = (40, 80, 120)

  mesh = RectilinearGrid3D((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :meg6)
  domain = mesh.iterators.cell.domain

  cell_volume = 0.05 * 0.025 * 0.025

  @test all(mesh.cell_center_metrics.J[domain] .≈ cell_volume)

  @test all(mesh.cell_center_metrics.x₁.ξ[domain] .≈ 0.05)
  @test all(mesh.cell_center_metrics.x₂.ξ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₃.ξ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₁.η[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₂.η[domain] .≈ 0.025)
  @test all(mesh.cell_center_metrics.x₃.η[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₁.ζ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₂.ζ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.x₃.ζ[domain] .≈ 0.025)

  @test all(mesh.cell_center_metrics.ξ̂.x₁[domain] .≈ 0.000625)
  @test all(mesh.cell_center_metrics.ξ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ξ̂.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.η̂.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.η̂.x₂[domain] .≈ 0.00125)
  @test all(mesh.cell_center_metrics.η̂.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ̂.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ̂.x₃[domain] .≈ 0.00125)

  @test all(mesh.cell_center_metrics.ξ.x₁[domain] .≈ 20.0)
  @test all(mesh.cell_center_metrics.ξ.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ξ.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.x₂[domain] .≈ 40.0)
  @test all(mesh.cell_center_metrics.η.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₃[domain] .≈ 40.0)
  @test all(mesh.cell_center_metrics.ξ.t[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.t[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.ζ.t[domain] .≈ 0.0)

  iaxis, jaxis, kaxis = (1, 2, 3)
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

  gcl(mesh)
end

@testset "3D Wavy Mesh GCL" begin
  using CurvilinearGrids
  using WriteVTK

  function wavy_grid(ni, nj, nk)
    Lx = Ly = Lz = 4.0

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    x = zeros(ni, nj, nk)
    y = zeros(ni, nj, nk)
    z = zeros(ni, nj, nk)
    for k in 1:nk
      for j in 1:nj
        for i in 1:ni
          x[i, j, k] = xmin + Δx0 * ((i - 1) + sinpi((j - 1) * Δy0) * sinpi((k - 1) * Δz0))
          y[i, j, k] = ymin + Δy0 * ((j - 1) + sinpi((k - 1) * Δz0) * sinpi((i - 1) * Δx0))
          z[i, j, k] = zmin + Δz0 * ((k - 1) + sinpi((i - 1) * Δx0) * sinpi((j - 1) * Δy0))
        end
      end
    end

    return (x, y, z)
  end

  ξ = 1
  η = 2
  ζ = 3

  ni = nj = nk = 20
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = CurvilinearGrid3D(x, y, z, :meg6)

  surf = extract_surface_mesh(mesh, :khi)
  save_vtk(mesh, "wavy3d")
  save_vtk(surf, "khi_surf")

  gcl(mesh)
end

@testset "3D Sphere Sector, Symmetric Conservative Metrics" begin
  r0, r1 = (1, 3)
  (θ0, θ1) = deg2rad.((35, 180 - 35))
  (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))

  ni, nj, nk = (20, 20, 20)
  mesh = rthetaphi_grid((r0, θ0, ϕ0), (r1, θ1, ϕ1), (ni, nj, nk), :meg6_symmetric)
  # mesh = rthetaphi_grid((r0, θ0, ϕ0), (r1, θ1, ϕ1), (ni, nj, nk), :meg6)
  save_vtk(mesh, "sphere_sector_3d")

  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  gcl(mesh)

  @test all(abs.(mesh.cell_center_metrics.ζ.x₃[mesh.iterators.cell.domain]) .< 5e-12)

  # @show I₁_max
  # @show I₂_max
  # @show I₃_max

  nothing
end
