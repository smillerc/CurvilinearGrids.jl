@testset "2D Rectangular Mesh" begin
  include("common.jl")

  ni, nj = (40, 80)
  nhalo = 5
  x0, x1 = (0, 2)
  y0, y1 = (1, 3)
  mesh = RectlinearGrid((x0, y0), (x1, y1), (ni, nj), nhalo)
  domain = mesh.iterators.cell.domain

  @test mesh.iterators.cell.full == CartesianIndices((50, 90))
  @test mesh.iterators.cell.domain == CartesianIndices((6:45, 6:85))

  @test mesh.iterators.node.full == CartesianIndices((51, 91))
  @test mesh.iterators.node.domain == CartesianIndices((6:46, 6:86))

  @test mesh.domain_limits.node == (ilo=6, ihi=46, jlo=6, jhi=86)
  @test mesh.domain_limits.cell == (ilo=6, ihi=45, jlo=6, jhi=85)

  @test all(mesh.cell_center_metrics.forward.J[domain] .≈ 0.00125)
  @test all(mesh.cell_center_metrics.forward.x₁.ξ[domain] .≈ 0.05)
  @test all(mesh.cell_center_metrics.forward.x₂.ξ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₁.η[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₂.η[domain] .≈ 0.025)

  @test all(mesh.cell_center_metrics.inverse.ξ.x₁[domain] .≈ 20.0)
  @test all(mesh.cell_center_metrics.inverse.ξ.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.η.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.η.x₂[domain] .≈ 40.0)
  @test all(mesh.cell_center_metrics.inverse.ξ.t[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.η.t[domain] .≈ 0.0)

  iaxis, jaxis = (1, 2)
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)

  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₁[i₊½_domain] .≈ 0.025)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₂[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.η̂.x₁[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.η̂.x₂[i₊½_domain] .≈ 0.05)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ξ̂.x₁[j₊½_domain] .≈ 0.025)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ξ̂.x₂[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.η̂.x₁[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.η̂.x₂[j₊½_domain] .≈ 0.05)

  @test all(mesh.edge_metrics.inverse.i₊½.ξ.x₁[i₊½_domain] .≈ 20.0)
  @test all(mesh.edge_metrics.inverse.i₊½.ξ.x₂[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse.i₊½.η.x₁[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse.i₊½.η.x₂[i₊½_domain] .≈ 40.0)
  @test all(mesh.edge_metrics.inverse.j₊½.ξ.x₁[j₊½_domain] .≈ 20.0)
  @test all(mesh.edge_metrics.inverse.j₊½.ξ.x₂[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse.j₊½.η.x₁[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse.j₊½.η.x₂[j₊½_domain] .≈ 40.0)

  ilo_c = nhalo + 1
  jlo_c = nhalo + 1

  @test coord(mesh, (ilo_c + 1, jlo_c + 1)) == [0.05, 1.025]
  @test centroid(mesh, (ilo_c, jlo_c)) == [0.025, 1.0125]

  xn, yn = coords(mesh)
  @test size(xn) == (41, 81)
  @test size(yn) == (41, 81)

  xc, yc = centroids(mesh)
  @test size(xc) == (40, 80)
  @test size(yc) == (40, 80)
end

@testset "2D Wavy Mesh GCL" begin
  include("common.jl")

  function wavy_grid(nx, ny)
    x0, x1 = (0, 1)
    y0, y1 = (0, 1)
    a0 = 0.1

    function x(i, j)
      x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
      y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
      return x1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
    end

    function y(i, j)
      x1d = x0 + (x1 - x0) * ((i - 1) / (nx - 1))
      y1d = y0 + (y1 - y0) * ((j - 1) / (ny - 1))
      return y1d + a0 * sin(2 * pi * x1d) * sin(2 * pi * y1d)
    end

    X = zeros(nx, ny)
    Y = zeros(nx, ny)
    for j in 1:ny
      for i in 1:nx
        X[i, j] = x(i, j)
        Y[i, j] = y(i, j)
      end
    end

    return X, Y
  end

  ni, nj = (41, 41)
  nhalo = 1
  x, y = wavy_grid(ni, nj)
  mesh = CurvilinearGrid2D(x, y, nhalo)

  domain = mesh.iterators.cell.domain

  ϵ = 5e-15
  I₁_passes = true
  I₂_passes = true
  for idx in domain
    i, j = idx.I

    ξ̂_i₊½ = mesh.edge_metrics.inverse_normalized.i₊½.ξ̂[i, j]
    ξ̂_i₋½ = mesh.edge_metrics.inverse_normalized.i₊½.ξ̂[i - 1, j]
    η̂_j₊½ = mesh.edge_metrics.inverse_normalized.j₊½.η̂[i, j]
    η̂_j₋½ = mesh.edge_metrics.inverse_normalized.j₊½.η̂[i, j - 1]

    I₁ = (ξ̂_i₊½.x₁ - ξ̂_i₋½.x₁) + (η̂_j₊½.x₁ - η̂_j₋½.x₁)
    I₂ = (ξ̂_i₊½.x₂ - ξ̂_i₋½.x₂) + (η̂_j₊½.x₂ - η̂_j₋½.x₂)

    I₁_passes = abs(I₁) < ϵ
    I₂_passes = abs(I₂) < ϵ
    if !(I₁_passes && I₂_passes)
      @show I₁ I₂
      break
    end
  end
  @test I₁_passes
  @test I₂_passes

  nothing

  save_vtk(mesh, "wavy")
end
