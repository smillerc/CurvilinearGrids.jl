
# @testitem "2D Mesh - Rectlinear Mesh" begin
include("common.jl")

function rect_grid(nx, ny)
  x0, x1 = (0, 2)
  y0, y1 = (1, 3)

  x(ξ, η) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
  y(ξ, η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
  return (x, y)
end

ni, nj = (5, 9)
nhalo = 0
x, y = rect_grid(ni, nj)
mesh = CurvilinearGrid2D(x, y, (ni, nj), nhalo)

plain_metrics = metrics(mesh, (2, 3))
@test metrics(mesh, (2, 3)) == (ξx=2.0, ξy=-0.0, ηx=-0.0, ηy=4.0, ξt=0.0, ηt=0.0)

@test conservative_metrics(mesh, (2, 3)) ==
  (ξ̂x=16.0, ξ̂y=-0.0, η̂x=-0.0, η̂y=32.0, ξt=0.0, ηt=0.0)

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

@benchmark conservative_metrics($mesh, (2, 3))
@benchmark inv(jacobian_matrix($mesh, (2, 3)))

bm1 = @benchmark metrics($mesh, (2, 3))
@test bm1.allocs == 0

bm2 = @benchmark jacobian_matrix($mesh, (2, 2))
@test bm2.allocs == 0

bm3 = @benchmark jacobian($mesh, (2, 2))
@test bm3.allocs == 0

ilo, ihi, jlo, jhi = mesh.limits

@test ilo = jlo == 1
@test ihi == 4
@test jhi == 8

@test coord(mesh, (1, 1)) == [0, 1]
@test coord(mesh, (2, 2)) == [0.5, 1.25]
@test coord(mesh, (2.5, 2.5)) == [0.75, 1.375]
@test coord(mesh, (2.5, 2.5)) == centroid(mesh, (2, 2))
@test centroid(mesh, (1, 1)) == [0.25, 1.125]

xy_coords = coords(mesh)
centroid_coords = centroids(mesh)
@test size(centroid_coords) == (2, 4, 8)
@test size(xy_coords) == (2, 5, 9)
# end
