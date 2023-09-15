
@testitem "2D Mesh - Rectlinear Mesh" begin
  include("common.jl")

  function rect_grid(nx, ny)
    x0, x1 = (0, 2)
    y0, y1 = (1, 3)

    _x(ξ) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
    _y(η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
    x(ξ, η) = _x(ξ)
    y(ξ, η) = _y(η)
    x((ξ, η)) = _x(ξ)
    y((ξ, η)) = _y(η)

    return (x, y)
  end

  ni, nj = (5, 9)
  nhalo = 0
  x, y = rect_grid(ni, nj)
  mesh = CurvilinearMesh2D(x, y, (ni, nj), nhalo)

  metrics(mesh, (2, 3)) == (ξ̂x=0.03125, ξ̂y=0.0, η̂x=0.0, η̂y=0.0625)
  bm1 = @benchmark metrics($mesh, (2, 3))
  @test bm1.allocs == 0

  @test jacobian_matrix(mesh, 2, 2) == @SMatrix [
    0.5 0.0
    0.0 0.25
  ]
  bm2 = @benchmark jacobian_matrix($mesh, (2, 2))
  @test bm2.allocs == 0

  cell_volume = 0.5 * 0.25
  @test jacobian(mesh, (2, 2)) == cell_volume

  bm3 = @benchmark jacobian($mesh, (2, 2))
  @test bm3.allocs == 0

  fn = "test2d"
  to_vtk(mesh, fn)
  @test isfile("$fn.vts")

  if isfile("$fn.vts")
    rm("$fn.vts")
  end

  ilo, ihi, jlo, jhi = mesh.limits

  @test ilo = jlo == 1
  @test ihi == 4
  @test jhi == 8

  for j in jlo:jhi
    for i in ilo:ihi
      m = metrics(mesh, (i, j))
      @test m == (ξ̂x=0.03125, ξ̂y=0.0, η̂x=0.0, η̂y=0.0625)
    end
  end

  @test coord(mesh, 1, 1) == [0, 1]
  @test coord(mesh, 2, 2) == [0.5, 1.25]
  @test coord(mesh, 2.5, 2.5) == [0.75, 1.375]
  @test coord(mesh, 2.5, 2.5) == centroid(mesh, 2, 2)
  @test centroid(mesh, 1, 1) == [0.25, 1.125]

  xy_coords = coords(mesh)
  centroid_coords = centroids(mesh)
  @test size(centroid_coords) == (2, 4, 8)
  @test size(xy_coords) == (2, 5, 9)
end
