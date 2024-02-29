
@testset "2D Rectangular mesh" begin
  function rect_grid(nx, ny)
    x0, x1 = (0, 2)
    y0, y1 = (1, 3)

    x(ξ, η) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
    y(ξ, η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
    return (x, y)
  end

  ni, nj = (5, 9)
  nhalo = 2
  x, y = rect_grid(ni, nj)
  mesh = CurvilinearGrid2D(x, y, (ni, nj), nhalo)
  domain = mesh.iterators.cell.domain

  @test mesh.iterators.cell.full == CartesianIndices((8, 12))
  @test mesh.iterators.cell.domain == CartesianIndices((3:6, 3:10))

  @test mesh.iterators.node.full == CartesianIndices((9, 13))
  @test mesh.iterators.node.domain == CartesianIndices((3:7, 3:11))

  @test mesh.domain_limits.node == (ilo=3, ihi=7, jlo=3, jhi=11)
  @test mesh.domain_limits.cell == (ilo=3, ihi=6, jlo=3, jhi=10)

  @test all(mesh.cell_center_metrics.J .≈ 0.125)
  @test all(mesh.cell_center_metrics.ξ.x .≈ 2.0)
  @test all(mesh.cell_center_metrics.ξ.y .≈ -0.0)
  @test all(mesh.cell_center_metrics.η.x .≈ -0.0)
  @test all(mesh.cell_center_metrics.η.y .≈ 4.0)
  @test all(mesh.cell_center_metrics.ξ.t .≈ 0.0)
  @test all(mesh.cell_center_metrics.η.t .≈ 0.0)

  @test all(mesh.edge_metrics.i₊½.J .≈ 0.125)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x .≈ 0.25)
  @test all(mesh.edge_metrics.i₊½.ξ̂.y .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.x .≈ 0.0)
  @test all(mesh.edge_metrics.i₊½.η̂.y .≈ 0.5)

  @test jacobian_matrix(mesh, (2, 2)) == @SMatrix [
    0.5 0.0
    0.0 0.25
  ]

  cell_area = 0.5 * 0.25
  @test jacobian(mesh, (2, 3)) == cell_area

  @test inv(jacobian_matrix(mesh, (2, 3))) == @SMatrix [
    2.0 0.0
    0.0 4.0
  ]

  @test inv(jacobian(mesh, (2, 3))) == 1 / cell_area

  bm0 = @benchmark conservative_metrics($mesh, (2, 3))
  @test bm0.allocs == 0

  bm1 = @benchmark metrics($mesh, (2, 3), 0)
  @test bm1.allocs == 0

  bm2 = @benchmark jacobian_matrix($mesh, (2, 2))
  @test bm2.allocs == 0

  bm3 = @benchmark jacobian($mesh, (2, 2))
  @test bm3.allocs == 0

  ilo_c = 3
  jlo_c = 3
  coord(mesh, (3, 3)) == [0, 1]
  @test coord(mesh, (ilo_c + 1, jlo_c + 1)) == [0.5, 1.25]
  @test coord(mesh, (ilo_c + 1.5, jlo_c + 1.5)) == [0.75, 1.375]
  @test coord(mesh, (2.5, 2.5)) == centroid(mesh, (2, 2))
  @test centroid(mesh, (ilo_c, jlo_c)) == [0.25, 1.125]

  xn, yn = coords(mesh)
  @test size(xn) == (5, 9)
  @test size(yn) == (5, 9)

  xc, yc = centroids(mesh)
  @test size(xc) == (4, 8)
  @test size(yc) == (4, 8)
end

@testset "2D Wavy Mesh GCL" begin
  # using Plots
  include("common.jl")
  # begin
  # include("../../src/grids/finitediff_metrics.jl")

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

    return (x, y)
  end

  ni, nj = (41, 41)
  nhalo = 1
  x, y = wavy_grid(ni, nj)
  mesh = CurvilinearGrid2D(x, y, (ni, nj), nhalo)

  domain = mesh.iterators.cell.domain

  I₁_passes = true
  I₂_passes = true
  for idx in domain
    i, j = idx.I .+ 0.5
    m_i₊½ = conservative_metrics(mesh, (i + 0.5, j))
    m_j₊½ = conservative_metrics(mesh, (i, j + 0.5))
    m_i₋½ = conservative_metrics(mesh, (i - 0.5, j))
    m_j₋½ = conservative_metrics(mesh, (i, j - 0.5))

    I₁ = (m_i₊½.ξ̂.x - m_i₋½.ξ̂.x) + (m_j₊½.η̂.x - m_j₋½.η̂.x)
    I₂ = (m_i₊½.ξ̂.y - m_i₋½.ξ̂.y) + (m_j₊½.η̂.y - m_j₋½.η̂.y)

    I₁_passes = abs(I₁) < eps()
    I₂_passes = abs(I₂) < eps()
    if !(I₁_passes && I₂_passes)
      break
    end
  end
  @test I₁_passes
  @test I₂_passes

  nothing

  # display(heatmap([ξ.x for ξ in mesh.cell_center_metrics.ξ]; title="ξx"))
  # display(heatmap([ξ.y for ξ in mesh.cell_center_metrics.ξ]; title="ξy"))
  # display(heatmap([ξ.x for ξ in mesh.cell_center_metrics.η]; title="ηx"))
  # display(heatmap([ξ.y for ξ in mesh.cell_center_metrics.η]; title="ηy"))
  # display(heatmap(mesh.cell_center_metrics.J; title="J"))
end
