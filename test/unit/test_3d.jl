using Test

@testset "3D Mesh - Rectangular Mesh" begin
  include("common.jl")

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
  nhalo = 5
  x, y, z = rect_grid(ni, nj, nk)

  mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
  domain = mesh.iterators.cell.domain

  cell_volume = 0.5 * 0.25 * 0.25

  @test all(mesh.cell_center_metrics.J[domain] .== cell_volume)
  @test all(mesh.cell_center_metrics.ξx[domain] .== 2.0)
  @test all(mesh.cell_center_metrics.ξy[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ξz[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ηx[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ηy[domain] .== 4.0)
  @test all(mesh.cell_center_metrics.ηz[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζx[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζy[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζz[domain] .== 4.0)
  @test all(mesh.cell_center_metrics.ξt[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ηt[domain] .== 0.0)
  @test all(mesh.cell_center_metrics.ζt[domain] .== 0.0)

  @show extrema(mesh.edge_metrics.i₊½.ξ̂x[domain])
  @test all(mesh.edge_metrics.i₊½.ξ̂x[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.j₊½.ξ̂x[domain] .≈ 0.0625)
  @test all(mesh.edge_metrics.k₊½.ξ̂x[domain] .≈ 0.0625)
  # @test all(mesh.edge_metrics.i₊½.ξ̂y[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.ξ̂y[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.ξ̂y[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.ξ̂z[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.ξ̂z[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.ξ̂z[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.η̂x[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.η̂x[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.η̂x[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.η̂y[domain] .≈ 0.125)
  # @test all(mesh.edge_metrics.j₊½.η̂y[domain] .≈ 0.125)
  # @test all(mesh.edge_metrics.k₊½.η̂y[domain] .≈ 0.125)
  # @test all(mesh.edge_metrics.i₊½.η̂z[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.η̂z[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.η̂z[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.ζ̂x[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.ζ̂x[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.ζ̂x[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.ζ̂y[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.ζ̂y[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.ζ̂y[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.ζ̂z[domain] .≈ 0.125)
  # @test all(mesh.edge_metrics.j₊½.ζ̂z[domain] .≈ 0.125)
  # @test all(mesh.edge_metrics.k₊½.ζ̂z[domain] .≈ 0.125)
  # @test all(mesh.edge_metrics.i₊½.ξ̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.ξ̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.ξ̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.η̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.η̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.η̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.i₊½.ζ̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.j₊½.ζ̂t[domain] .≈ 0.0)
  # @test all(mesh.edge_metrics.k₊½.ζ̂t[domain] .≈ 0.0)

  conserved_metrics_pass = false
  for idx in domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂x[i, j, k] - mesh.edge_metrics.i₊½.ξ̂x[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂x[i, j, k] - mesh.edge_metrics.j₊½.η̂x[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂x[i, j, k] - mesh.edge_metrics.k₊½.ζ̂x[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂y[i, j, k] - mesh.edge_metrics.i₊½.ξ̂y[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂y[i, j, k] - mesh.edge_metrics.j₊½.η̂y[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂y[i, j, k] - mesh.edge_metrics.k₊½.ζ̂y[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂z[i, j, k] - mesh.edge_metrics.i₊½.ξ̂z[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂z[i, j, k] - mesh.edge_metrics.j₊½.η̂z[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂z[i, j, k] - mesh.edge_metrics.k₊½.ζ̂z[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > 4eps())
    I₂ = I₂ * (abs(I₂) > 4eps())
    I₃ = I₃ * (abs(I₃) > 4eps())

    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  # @test conserved_metrics_pass

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

@testset "3D Wavy Mesh Conserved Metrics GCL" begin
  using CurvilinearGrids
  using CurvilinearGrids.MetricDiscretizationSchemes.MonotoneExplicit6thOrderScheme:
    conserved_metric!

  include("../../src/metric_schemes/indexing_fun.jl")

  @inline function ∂(ϕ, idx, dim)
    @inbounds begin
      f1 = 3 / 4
      f2 = 3 / 20
      f3 = 1 / 60

      ᵢ₊₁ = up(idx, dim, 1)
      ᵢ₊₂ = up(idx, dim, 2)
      ᵢ₊₃ = up(idx, dim, 3)
      ᵢ₋₁ = down(idx, dim, 1)
      ᵢ₋₂ = down(idx, dim, 2)
      ᵢ₋₃ = down(idx, dim, 3)

      ϕᵢ₊₁ = ϕ(ᵢ₊₁)
      ϕᵢ₊₂ = ϕ(ᵢ₊₂)
      ϕᵢ₊₃ = ϕ(ᵢ₊₃)
      ϕᵢ₋₁ = ϕ(ᵢ₋₁)
      ϕᵢ₋₂ = ϕ(ᵢ₋₂)
      ϕᵢ₋₃ = ϕ(ᵢ₋₃)

      ∂ϕ = (f1 * (ϕᵢ₊₁ - ϕᵢ₋₁) - f2 * (ϕᵢ₊₂ - ϕᵢ₋₂) + f3 * (ϕᵢ₊₃ - ϕᵢ₋₃))
    end
    return ∂ϕ
  end

  @inline function ∂²(ϕ, idx, dim)
    @inbounds begin
      f1 = 1 / 2

      ᵢ₊₁ = up(idx, dim, 1)
      ᵢ₋₁ = down(idx, dim, 1)
      ϕᵢ₊₁ = ϕ(ᵢ₊₁)
      ϕᵢ = ϕ(idx)
      ϕᵢ₋₁ = ϕ(ᵢ₋₁)

      ∂ϕᵢ₊₁ = ∂(ϕ, ᵢ₊₁, dim)
      ∂ϕᵢ₋₁ = ∂(ϕ, ᵢ₋₁, dim)

      ∂²ϕ = (2(ϕᵢ₊₁ - 2ϕᵢ + ϕᵢ₋₁) - f1 * (∂ϕᵢ₊₁ - ∂ϕᵢ₋₁))
    end
    return ∂²ϕ
  end

  function ∂ϕ(ϕ, idx, dim)
    @inbounds begin
      ᵢ₊₁ = up(idx, dim, 1)
      ᵢ₋₁ = down(idx, dim, 1)

      xᵢ = ϕ(idx)
      xᵢ₊₁ = ϕ(ᵢ₊₁)
      xᵢ₋₁ = ϕ(ᵢ₋₁)

      xᴸᵢ₊½ = xᵢ + 0.5∂(ϕ, idx, dim) + ∂²(ϕ, idx, dim) / 12
      xᴿᵢ₊½ = xᵢ₊₁ - 0.5∂(ϕ, ᵢ₊₁, dim) + ∂²(ϕ, ᵢ₊₁, dim) / 12
      xᵢ₊½ = (xᴿᵢ₊½ + xᴸᵢ₊½) / 2

      xᴸᵢ₋½ = xᵢ₋₁ + 0.5∂(ϕ, ᵢ₋₁, dim) + ∂²(ϕ, ᵢ₋₁, dim) / 12
      xᴿᵢ₋½ = xᵢ - 0.5∂(ϕ, idx, dim) + ∂²(ϕ, idx, dim) / 12

      xᵢ₋½ = (xᴿᵢ₋½ + xᴸᵢ₋½) / 2
    end
    return xᵢ₊½ - xᵢ₋½
  end

  function ϕᵢ₊½(ϕ, idx, dim)
    @inbounds begin
      ᵢ₊₁ = up(idx, dim, 1)
      xᵢ = ϕ(idx)
      xᵢ₊₁ = ϕ(ᵢ₊₁)
      xᴸᵢ₊½ = xᵢ + 0.5∂(ϕ, idx, dim) + ∂²(ϕ, idx, dim) / 12
      xᴿᵢ₊½ = xᵢ₊₁ - 0.5∂(ϕ, ᵢ₊₁, dim) + ∂²(ϕ, ᵢ₊₁, dim) / 12
      xᵢ₊½ = (xᴿᵢ₊½ + xᴸᵢ₊½) / 2
    end
    return xᵢ₊½
  end

  function save_vtk(mesh)
    fn = "wavy"
    @info "Writing to $fn.vti"

    xyz_n = CurvilinearGrids.coords(mesh)
    domain = mesh.iterators.cell.domain

    # @views begin
    #   #   J = [m.J for m in mesh.cell_center_metrics[domain]]
    #   ξx = [m.x for m in mesh.cell_center_metrics.ξ[domain]]
    #   ξ̂xᵢ₊½ = [m.x for m in mesh.edge_metrics.i₊½.ξ̂[domain]]
    #   #   ηy = [m.ηy for m in mesh.cell_center_metrics[domain]]
    # end

    vtk_grid(fn, xyz_n) do vtk

      # w1 = [x[1] for idx in eachindex(f.char_state.W[ilo:ihi, jlo:jhi]
      for (k, v) in pairs(mesh.cell_center_metrics)
        vtk["$k"] = v[domain]
      end

      for (edge_name, data) in pairs(mesh.edge_metrics)
        for (k, v) in pairs(data)
          vtk["$(k)_$(edge_name)"] = v[domain]
        end
      end
      # vtk["ξx"] = ξx
      # vtk["ξ̂xᵢ₊½"] = ξ̂xᵢ₊½
      # vtk["I"] = I[domain]
      # vtk["ηy"] = ηy
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
  nhalo = 5
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)

  xc(idx) = mesh.coord_funcs.x((idx.I .+ 0.5 .- nhalo)...)
  yc(idx) = mesh.coord_funcs.y((idx.I .+ 0.5 .- nhalo)...)
  zc(idx) = mesh.coord_funcs.z((idx.I .+ 0.5 .- nhalo)...)
  # xc(idx) = mesh.centroids.x[idx]
  # yc(idx) = mesh.centroids.y[idx]
  # zc(idx) = mesh.centroids.z[idx]

  xζ(idx) = ∂ϕ(xc, idx, ζ)
  xη(idx) = ∂ϕ(xc, idx, η)
  xξ(idx) = ∂ϕ(xc, idx, ξ)
  yζ(idx) = ∂ϕ(yc, idx, ζ)
  yη(idx) = ∂ϕ(yc, idx, η)
  yξ(idx) = ∂ϕ(yc, idx, ξ)
  zζ(idx) = ∂ϕ(zc, idx, ζ)
  zη(idx) = ∂ϕ(zc, idx, η)
  zξ(idx) = ∂ϕ(zc, idx, ξ)

  xζy(idx) = ∂ϕ(xc, idx, ζ) * yc(idx)
  xηy(idx) = ∂ϕ(xc, idx, η) * yc(idx)
  xξy(idx) = ∂ϕ(xc, idx, ξ) * yc(idx)
  yζz(idx) = ∂ϕ(yc, idx, ζ) * zc(idx)
  yηz(idx) = ∂ϕ(yc, idx, η) * zc(idx)
  yξz(idx) = ∂ϕ(yc, idx, ξ) * zc(idx)
  zζx(idx) = ∂ϕ(zc, idx, ζ) * xc(idx)
  zηx(idx) = ∂ϕ(zc, idx, η) * xc(idx)
  zξx(idx) = ∂ϕ(zc, idx, ξ) * xc(idx)

  xζy_η(idx) = ∂ϕ(xζy, idx, η)
  xζy_ξ(idx) = ∂ϕ(xζy, idx, ξ)
  xηy_ζ(idx) = ∂ϕ(xηy, idx, ζ)
  xηy_ξ(idx) = ∂ϕ(xηy, idx, ξ)
  xξy_ζ(idx) = ∂ϕ(xξy, idx, ζ)
  xξy_η(idx) = ∂ϕ(xξy, idx, η)
  yζz_η(idx) = ∂ϕ(yζz, idx, η)
  yζz_ξ(idx) = ∂ϕ(yζz, idx, ξ)
  yηz_ζ(idx) = ∂ϕ(yηz, idx, ζ)
  yηz_ξ(idx) = ∂ϕ(yηz, idx, ξ)
  yξz_ζ(idx) = ∂ϕ(yξz, idx, ζ)
  yξz_η(idx) = ∂ϕ(yξz, idx, η)
  zζx_η(idx) = ∂ϕ(zζx, idx, η)
  zζx_ξ(idx) = ∂ϕ(zζx, idx, ξ)
  zηx_ζ(idx) = ∂ϕ(zηx, idx, ζ)
  zηx_ξ(idx) = ∂ϕ(zηx, idx, ξ)
  zξx_ζ(idx) = ∂ϕ(zξx, idx, ζ)
  zξx_η(idx) = ∂ϕ(zξx, idx, η)

  ξ̂x(idx) = yηz_ζ(idx) - yζz_η(idx)
  ξ̂y(idx) = zηx_ζ(idx) - zζx_η(idx)
  ξ̂z(idx) = xηy_ζ(idx) - xζy_η(idx)

  η̂x(idx) = yζz_ξ(idx) - yξz_ζ(idx)
  η̂y(idx) = zζx_ξ(idx) - zξx_ζ(idx)
  η̂z(idx) = xζy_ξ(idx) - xξy_ζ(idx)

  ζ̂x(idx) = yξz_η(idx) - yηz_ξ(idx)
  ζ̂y(idx) = zξx_η(idx) - zηx_ξ(idx)
  ζ̂z(idx) = xξy_η(idx) - xηy_ξ(idx)

  ξ̂xᵢ₊½(idx) = ϕᵢ₊½(ξ̂x, idx, ξ)
  η̂xᵢ₊½(idx) = ϕᵢ₊½(η̂x, idx, ξ)
  ζ̂xᵢ₊½(idx) = ϕᵢ₊½(ζ̂x, idx, ξ)

  ξ̂yⱼ₊½(idx) = ϕᵢ₊½(ξ̂y, idx, η)
  η̂yⱼ₊½(idx) = ϕᵢ₊½(η̂y, idx, η)
  ζ̂yⱼ₊½(idx) = ϕᵢ₊½(ζ̂y, idx, η)

  ξ̂zₖ₊½(idx) = ϕᵢ₊½(ξ̂z, idx, ζ)
  η̂zₖ₊½(idx) = ϕᵢ₊½(η̂z, idx, ζ)
  ζ̂zₖ₊½(idx) = ϕᵢ₊½(ζ̂z, idx, ζ)

  ∂ξ̂x∂ξ(idx) = ∂ϕ(ξ̂x, idx, ξ)
  ∂η̂x∂η(idx) = ∂ϕ(η̂x, idx, η)
  ∂ζ̂x∂ζ(idx) = ∂ϕ(ζ̂x, idx, ζ)
  ∂ξ̂y∂ξ(idx) = ∂ϕ(ξ̂y, idx, ξ)
  ∂η̂y∂η(idx) = ∂ϕ(η̂y, idx, η)
  ∂ζ̂y∂ζ(idx) = ∂ϕ(ζ̂y, idx, ζ)
  ∂ξ̂z∂ξ(idx) = ∂ϕ(ξ̂z, idx, ξ)
  ∂η̂z∂η(idx) = ∂ϕ(η̂z, idx, η)
  ∂ζ̂z∂ζ(idx) = ∂ϕ(ζ̂z, idx, ζ)

  # save_vtk(mesh)
  domain = mesh.iterators.cell.full

  # i, j, k = (6, 16, 15)
  # idx = CartesianIndex(i, j, k)

  # ξ̂_x = similar(mesh.cell_center_metrics.J)
  # η_axis = 2
  # ζ_axis = 3

  # println("woohoo!")

  # conserved_metric!(
  #   mesh.discretization_scheme,
  #   ξ̂_x,
  #   mesh.centroids.y,
  #   η_axis,
  #   mesh.centroids.z,
  #   ζ_axis,
  #   domain,
  # )

  # _yη = @views mesh.discretization_scheme.cache.inner_deriv1
  # _yηz_ζ = @views mesh.discretization_scheme.cache.outer_deriv1
  # _yζ = @views mesh.discretization_scheme.cache.inner_deriv2
  # _yζz_η = @views mesh.discretization_scheme.cache.outer_deriv2

  # # @show domain
  # # for idx in domain
  # @show yη(idx) - _yη[idx]
  # @show yηz(idx) - _yη[idx] * zc(idx)
  # @show yηz_ζ(idx) - _yηz_ζ[idx]

  # # println()

  # d = yζ(idx) - _yζ[idx]
  # # @show yζz(idx) - _yζ[idx] * zc(idx)
  # # @show yζz_η(idx) - _yζz_η[idx] # here's where the error shows up for some reason!

  # # d = ξ̂x(idx) - ξ̂_x[idx]
  # if abs(d) > 1e-14
  #   @show idx, d
  #   # break
  # end
  # # end

  # # @show ξ̂xᵢ₊½(idx) - mesh.edge_metrics.i₊½.ξ̂x[idx]
  # # @show η̂xᵢ₊½(idx), mesh.edge_metrics.i₊½.η̂x[idx]

  conserved_metrics_pass = false
  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (mesh.edge_metrics.i₊½.ξ̂x[i, j, k] - mesh.edge_metrics.i₊½.ξ̂x[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂x[i, j, k] - mesh.edge_metrics.j₊½.η̂x[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂x[i, j, k] - mesh.edge_metrics.k₊½.ζ̂x[i, j, k - 1])
    )
    I₂ = (
      (mesh.edge_metrics.i₊½.ξ̂y[i, j, k] - mesh.edge_metrics.i₊½.ξ̂y[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂y[i, j, k] - mesh.edge_metrics.j₊½.η̂y[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂y[i, j, k] - mesh.edge_metrics.k₊½.ζ̂y[i, j, k - 1])
    )
    I₃ = (
      (mesh.edge_metrics.i₊½.ξ̂z[i, j, k] - mesh.edge_metrics.i₊½.ξ̂z[i - 1, j, k]) +
      (mesh.edge_metrics.j₊½.η̂z[i, j, k] - mesh.edge_metrics.j₊½.η̂z[i, j - 1, k]) +
      (mesh.edge_metrics.k₊½.ζ̂z[i, j, k] - mesh.edge_metrics.k₊½.ζ̂z[i, j, k - 1])
    )

    I₁ = I₁ * (abs(I₁) > 4eps())
    I₂ = I₂ * (abs(I₂) > 4eps())
    I₃ = I₃ * (abs(I₃) > 4eps())

    # @show I₁, I₂, I₃
    conserved_metrics_pass = iszero(I₁) && iszero(I₂) && iszero(I₃)
    if !conserved_metrics_pass
      break
    end
  end
  @test conserved_metrics_pass
end

@testset "3D Mesh - Sphere Sector" begin
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
  nhalo = 5
  x, y, z = sphere_grid(ni, nj, nk)
  @test_nowarn CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
end

nothing