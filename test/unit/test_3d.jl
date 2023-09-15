
# @testitem "3D Mesh - Rectlinear Mesh" begin
#   include("common.jl")

function rect_grid(nx, ny, nz)
  x0, x1 = (0, 2)
  y0, y1 = (1, 3)
  z0, z1 = (-1, 2)

  _x(ξ) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
  _y(η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
  _z(ζ) = @. z0 + (z1 - z0) * ((ζ - 1) / (nz - 1))

  x(ξ, η, ζ) = _x(ξ)
  y(ξ, η, ζ) = _y(η)
  z(ξ, η, ζ) = _z(ζ)

  x((ξ, η, ζ)) = _x(ξ)
  y((ξ, η, ζ)) = _y(η)
  z((ξ, η, ζ)) = _z(ζ)

  return (x, y, z)
end

ni, nj, nk = (5, 9, 8)
nhalo = 0
x, y, z = rect_grid(ni, nj, nk)
mesh = CurvilinearMesh3D(x, y, z, (ni, nj, nk), nhalo)

metrics(mesh, 2, 3, 4) # (2, 3, 4)
# to_vtk(mesh, "test3d")

# ilo, ihi, jlo, jhi, klo, khi = mesh.limits

# @test ilo = jlo == 1
# @test ihi == 4
# @test jhi == 8
# @test klo == 1
# @test khi == 7
# metric_tuple = mesh.metrics

# for j in jlo:jhi
#   for i in ilo:ihi
#     @test metric_tuple.∂x∂ξ(i, j) == 0.5
#     @test metric_tuple.∂x∂η(i, j) == 0.0
#     @test metric_tuple.∂y∂ξ(i, j) == 0.0
#     @test metric_tuple.∂y∂η(i, j) == 0.25
#   end
# end

# @test xy(mesh, 1, 1) == [0, 1]
# @test xy(mesh, 2, 2) == [0.5, 1.25]
# @test xy(mesh, 2.5, 2.5) == [0.75, 1.375]
# @test xy(mesh, 2.5, 2.5) == centroid_xy(mesh, 2, 2)
# @test centroid_xy(mesh, 1, 1) == [0.25, 1.125]

# xy_coords = coords(mesh)
# centroid_coords = centroids(mesh)
# @test size(centroid_coords) == (2, 4, 8)
# @test size(xy_coords) == (2, 5, 9)

# all([metric_tuple.J.(i, j) for i in ilo:ihi for j in jlo:jhi] .== 0.125)
# all([metric_tuple.J⁻¹.(i, j) for i in ilo:ihi for j in jlo:jhi] .== 8.0)
# end
