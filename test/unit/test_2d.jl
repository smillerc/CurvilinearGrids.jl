
@testitem "2D Mesh - Rectlinear Mesh" begin
  include("common.jl")

  function rect_grid(nx, ny)
    x0, x1 = (0, 2)
    y0, y1 = (1, 3)

    _x(ξ) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
    _y(η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
    x(ξ, η) = _x(ξ)
    y(ξ, η) = _y(η)

    return (x, y)
  end

  ni, nj = (5, 9)
  nhalo = 0
  x, y = rect_grid(ni, nj)
  mesh = CurvilinearMesh2D(x, y, (ni, nj), nhalo)

  ilo, ihi, jlo, jhi = mesh.limits

  @test ilo = jlo == 1
  @test ihi == 4
  @test jhi == 8
  metric_tuple = mesh.metrics

  for j in jlo:jhi
    for i in ilo:ihi
      @test metric_tuple.∂x∂ξ(i, j) == 0.5
      @test metric_tuple.∂x∂η(i, j) == 0.0
      @test metric_tuple.∂y∂ξ(i, j) == 0.0
      @test metric_tuple.∂y∂η(i, j) == 0.25
    end
  end

  @test xy(mesh, 1, 1) == [0, 1]
  @test xy(mesh, 2, 2) == [0.5, 1.25]
  @test xy(mesh, 2.5, 2.5) == [0.75, 1.375]
  @test xy(mesh, 2.5, 2.5) == centroid_xy(mesh, 2, 2)
  @test centroid_xy(mesh, 1, 1) == [0.25, 1.125]

  xy_coords = coords(mesh)
  centroid_coords = centroids(mesh)
  @test size(centroid_coords) == (2, 4, 8)
  @test size(xy_coords) == (2, 5, 9)

  all([metric_tuple.J.(i, j) for i in ilo:ihi for j in jlo:jhi] .== 0.125)
  all([metric_tuple.J⁻¹.(i, j) for i in ilo:ihi for j in jlo:jhi] .== 8.0)
end
