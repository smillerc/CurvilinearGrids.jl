
@testset "3D Rectangular Mesh Metrics, Conserved Metrics, GCL" begin
  function rect_grid(nx, ny, nz)
    x0, x1 = (0.0, 2.0)
    y0, y1 = (1, 3)
    z0, z1 = (-1, 2)

    x(ξ, η, ζ) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
    y(ξ, η, ζ) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
    z(ξ, η, ζ) = @. z0 + (z1 - z0) * ((ζ - 1) / (nz - 1))

    return (x, y, z)
  end

  ni, nj, nk = (5, 9, 13)
  nhalo = 6
  x, y, z = rect_grid(ni, nj, nk)

  mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
  domain = mesh.iterators.cell.domain

  cell_volume = 0.5 * 0.25 * 0.25

  @test all(mesh.cell_center_metrics.J[domain] .== cell_volume)
  @test all(mesh.cell_center_metrics.ξ.x[domain] .== 2.0)
  @test all(mesh.cell_center_metrics.ξ.y[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ξ.z[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.η.x[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.η.y[domain] .== 4.0)
  @test all(mesh.cell_center_metrics.η.z[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.x[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.y[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.z[domain] .== 4.0)
  @test all(mesh.cell_center_metrics.ξ.t[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.η.t[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.t[domain] .== 0.0)

  @test all(mesh.edge_metrics.i₊½.ξ̂.x[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.k₊½.ξ̂.x[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.i₊½.ξ̂.y[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.y[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.y[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ξ̂.z[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.z[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.z[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.x[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.x[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.x[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.y[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.j₊½.η̂.y[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.k₊½.η̂.y[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.i₊½.η̂.z[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.z[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.z[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.x[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.x[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.x[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.y[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.y[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.y[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.z[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.j₊½.ζ̂.z[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.k₊½.ζ̂.z[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.i₊½.ξ̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.t[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.t[domain] .≈ 0.0)

  conserved_metrics_pass = false
  # for idx in expand(domain, -1)
  for idx in domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂.x[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x[i, j, k] - mesh.edge_metrics.j₊½.η̂.x[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂.y[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.y[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.y[i, j, k] - mesh.edge_metrics.j₊½.η̂.y[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.y[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.y[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂.z[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.z[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.z[i, j, k] - mesh.edge_metrics.j₊½.η̂.z[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.z[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.z[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > 4eps())
    I₂ = I₂ * (abs(I₂) > 4eps())
    I₃ = I₃ * (abs(I₃) > 4eps())

    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  @test conserved_metrics_pass

  # bm1 = @benchmark metrics($mesh, (2, 3, 4))
  # @test bm1.allocs == 0

  # @test jacobian_matrix(mesh, (2, 3, 4)) == @SMatrix [
  #   0.5 0.0 0.0
  #   0.0 0.25 0.0
  #   0.0 0.0 0.25
  # ]

  # bm2 = @benchmark jacobian_matrix($mesh, (2, 3, 4.5))
  # @test bm2.allocs == 0

  # cell_volume = 0.5 * 0.25 * 0.25
  # @test jacobian(mesh, (2, 3, 4)) == cell_volume

  # bm3 = @benchmark jacobian($mesh, (2, 3, 4))
  # @test bm3.allocs == 0

  # ilo, ihi, jlo, jhi, klo, khi = mesh.limits
  # @test ilo = jlo = klo == 1
  # @test ihi == 4
  # @test jhi == 8
  # @test khi == 12

  # # for k in klo:khi
  # #   for j in jlo:jhi
  # #     for i in ilo:ihi
  # #       @test metrics(mesh, (i, j, k)) == m
  # #     end
  # #   end
  # # end

  # @test coord(mesh, (1, 1, 1)) == [0, 1, -1]
  # @test coord(mesh, (2.5, 2.5, 2.5)) == [0.75, 1.375, -0.625]
  # @test coord(mesh, (2.5, 2.5, 2.5)) == centroid(mesh, (2, 2, 2))

  # xyz_coords = coords(mesh)
  # centroid_coords = centroids(mesh)
  # @test size(centroid_coords) == (3, 4, 8, 12)
  # @test size(xyz_coords) == (3, 5, 9, 13)
end

@testset "3D Wavy Mesh GCL" begin
  using CurvilinearGrids
  using WriteVTK
  using CurvilinearGrids.MetricDiscretizationSchemes.MonotoneExplicit6thOrderScheme:
    conserved_metric!

  # include("../../src/metric_schemes/indexing_fun.jl")

  function save_vtk(mesh)
    fn = "wavy"
    @info "Writing to $fn.vti"

    xyz_n = CurvilinearGrids.coords(mesh)
    domain = mesh.iterators.cell.domain

    @views vtk_grid(fn, xyz_n) do vtk
      for (k, v) in pairs(mesh.cell_center_metrics)
        vtk["$k"] = v[domain]
      end

      for (edge_name, data) in pairs(mesh.edge_metrics)
        for (k, v) in pairs(data)
          vtk["$(k)_$(edge_name)"] = v[domain]
        end
      end
    end
  end

  function wavy_grid(ni, nj, nk)
    Lx = Ly = Lz = 4.0

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    x(i, j, k) = xmin + Δx0 * ((i - 1) + sinpi((j - 1) * Δy0) * sinpi((k - 1) * Δz0))
    y(i, j, k) = ymin + Δy0 * ((j - 1) + sinpi((k - 1) * Δz0) * sinpi((i - 1) * Δx0))
    z(i, j, k) = zmin + Δz0 * ((k - 1) + sinpi((i - 1) * Δx0) * sinpi((j - 1) * Δy0))

    return (x, y, z)
  end

  ξ = 1
  η = 2
  ζ = 3

  ni = nj = nk = 20
  nhalo = 6
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)

  # save_vtk(mesh)

  conserved_metrics_pass = false
  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂.x[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x[i, j, k] - mesh.edge_metrics.j₊½.η̂.x[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂.y[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.y[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.y[i, j, k] - mesh.edge_metrics.j₊½.η̂.y[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.y[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.y[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂.z[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.z[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.z[i, j, k] - mesh.edge_metrics.j₊½.η̂.z[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.z[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.z[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > 10eps())
    I₂ = I₂ * (abs(I₂) > 10eps())
    I₃ = I₃ * (abs(I₃) > 10eps())

    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  @test conserved_metrics_pass
end

@testset "3D Sphere Sector Mesh Construction" begin
  include("common.jl")

  function sphere_grid(nr, ntheta, nphi)
    r0, r1 = (1, 3)
    (θ0, θ1) = deg2rad.((35, 180 - 35))
    (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))

    r(ξ) = r0 + (r1 - r0) * ((ξ - 1) / (nr - 1))
    θ(η) = θ0 + (θ1 - θ0) * ((η - 1) / (ntheta - 1))
    ϕ(ζ) = ϕ0 + (ϕ1 - ϕ0) * ((ζ - 1) / (nphi - 1))

    x(ξ, η, ζ) = r(ξ) * sin(θ(η)) * cos(ϕ(ζ))
    y(ξ, η, ζ) = r(ξ) * sin(θ(η)) * sin(ϕ(ζ))
    z(ξ, η, ζ) = r(ξ) * cos(θ(η))

    return (x, y, z)
  end

  ni, nj, nk = (5, 9, 11)
  nhalo = 6
  x, y, z = sphere_grid(ni, nj, nk)
  @test_nowarn CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
end
