
@testset "3D Rectangular Mesh Metrics, Conserved Metrics, GCL" begin
  include("common.jl")

  x0, x1 = (0.0, 2.0)
  y0, y1 = (1, 3)
  z0, z1 = (-1, 2)
  ni, nj, nk = (4, 8, 12)
  nhalo = 4

  mesh = RectlinearGrid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), nhalo)
  domain = mesh.iterators.cell.domain

  cell_volume = 0.5 * 0.25 * 0.25

  @test all(mesh.cell_center_metrics.J[domain] .== cell_volume)
  @test all(mesh.cell_center_metrics.ξ.x₁[domain] .== 2.0)
  @test all(mesh.cell_center_metrics.ξ.x₂[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ξ.x₃[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.η.x₁[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.η.x₂[domain] .== 4.0)
  @test all(mesh.cell_center_metrics.η.x₃[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₁[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₂[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.x₃[domain] .== 4.0)
  @test all(mesh.cell_center_metrics.ξ.t[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.η.t[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζ.t[domain] .== 0.0)

  @test all(mesh.edge_metrics.i₊½.ξ̂.x₁[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₁[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.k₊½.ξ̂.x₁[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₃[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ξ̂.x₃[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ξ̂.x₃[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.x₁[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.x₁[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.x₁[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.x₂[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.j₊½.η̂.x₂[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.k₊½.η̂.x₂[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.i₊½.η̂.x₃[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.η̂.x₃[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.η̂.x₃[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.x₁[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.x₁[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.x₁[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.j₊½.ζ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.k₊½.ζ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.ζ̂.x₃[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.j₊½.ζ̂.x₃[domain] .≈ 0.125)
  @test all(mesh.edge_metrics.k₊½.ζ̂.x₃[domain] .≈ 0.125)
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
      (mesh.edge_metrics.i₊½.ξ̂.x₁[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₁[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₁[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₂[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₂[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₂[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₃[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₃[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₃[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k - 1])
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
  # using CurvilinearGrids.MetricDiscretizationSchemes.MonotoneExplicit6thOrderScheme:
  #   conserved_metric!

  # include("../../src/metric_schemes/indexing_fun.jl")

  # function save_vtk(mesh)
  #   fn = "wavy"
  #   @info "Writing to $fn.vti"

  #   xyz_n = CurvilinearGrids.coords(mesh)
  #   domain = mesh.iterators.cell.domain

  #   @views vtk_grid(fn, xyz_n) do vtk
  #     for (k, v) in pairs(mesh.cell_center_metrics)
  #       vtk["$k"] = v[domain]
  #     end

  #     for (edge_name, data) in pairs(mesh.edge_metrics)
  #       for (k, v) in pairs(data)
  #         vtk["$(k)_$(edge_name)"] = v[domain]
  #       end
  #     end
  #   end
  # end

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
  nhalo = 4
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = CurvilinearGrid3D(x, y, z, nhalo)

  save_vtk(mesh, "wavy3d")

  conserved_metrics_pass = false
  ϵ = 1e-14
  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₁[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₁[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₁[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₂[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₂[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₂[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₃[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₃[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₃[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > ϵ)
    I₂ = I₂ * (abs(I₂) > ϵ)
    I₃ = I₃ * (abs(I₃) > ϵ)

    # @show I₁, I₂, I₃
    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  @test conserved_metrics_pass
end

@testset "3D Sphere Sector Mesh Construction" begin
  include("common.jl")

  r0, r1 = (1, 3)
  (θ0, θ1) = deg2rad.((35, 180 - 35))
  (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))
  ni, nj, nk = (20, 20, 20)
  nhalo = 4
  # @test_nowarn CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
  mesh = RThetaPhiGrid((r0, θ0, ϕ0), (r1, θ1, ϕ1), (ni, nj, nk), nhalo)
  save_vtk(mesh, "sphere_sector_3d")

  conserved_metrics_pass = false
  ϵ = 1e-14
  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₁[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₁[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₁[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₂[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₂[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₂[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂.x₃[i, j, k] - mesh.edge_metrics.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂.x₃[i, j, k] - mesh.edge_metrics.j₊½.η̂.x₃[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k] - mesh.edge_metrics.k₊½.ζ̂.x₃[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > ϵ)
    I₂ = I₂ * (abs(I₂) > ϵ)
    I₃ = I₃ * (abs(I₃) > ϵ)

    # @show I₁, I₂, I₃
    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  @test conserved_metrics_pass
end
