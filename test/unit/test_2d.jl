
# @testset "2D rectangular mesh" begin
include("common.jl")

function rect_grid(nx, ny)
  x0, x1 = (0, 2)
  y0, y1 = (1, 3)

  x(ξ, η) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))
  y(ξ, η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
  return (x, y)
end

ni, nj = (50, 90)
nhalo = 5
x, y = rect_grid(ni, nj)
mesh = CurvilinearGrid2D(x, y, (ni, nj), nhalo)
domain = mesh.iterators.cell.domain

cswh = cellsize_withhalo(mesh)
meg6_metrics = MEG6Scheme(cswh)

# centroids
xc = zeros(cswh)
yc = zeros(cswh)
for idx in CartesianIndices(xc)
  i, j = idx.I
  xc[idx] = mesh.x_func(i - mesh.nhalo + 0.5, j - mesh.nhalo + 0.5) # x
  yc[idx] = mesh.y_func(i - mesh.nhalo + 0.5, j - mesh.nhalo + 0.5) # y
end

domain = mesh.iterators.cell.domain

CurvilinearGrids.MetricDiscretizationSchemes.update_metrics!(meg6_metrics, xc, yc, domain)

meg6_metrics.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₁.i₊½

@test mesh.iterators.cell.full == CartesianIndices((8, 12))
@test mesh.iterators.cell.domain == CartesianIndices((3:6, 3:10))
@test mesh.iterators.cell.ilo_halo == CartesianIndices((1:2, 3:10))
@test mesh.iterators.cell.ihi_halo == CartesianIndices((7:8, 3:10))
@test mesh.iterators.cell.jlo_halo == CartesianIndices((3:6, 1:2))
@test mesh.iterators.cell.jhi_halo == CartesianIndices((3:6, 11:12))

@test mesh.iterators.node.full == CartesianIndices((9, 13))
@test mesh.iterators.node.domain == CartesianIndices((3:7, 3:11))
@test mesh.iterators.node.ilo_halo == CartesianIndices((1:3, 3:11))
@test mesh.iterators.node.ihi_halo == CartesianIndices((7:9, 3:11))
@test mesh.iterators.node.jlo_halo == CartesianIndices((3:7, 1:3))
@test mesh.iterators.node.jhi_halo == CartesianIndices((3:7, 11:13))

@test mesh.domain_limits.node == (ilo=3, ihi=7, jlo=3, jhi=11)
@test mesh.domain_limits.cell == (ilo=3, ihi=6, jlo=3, jhi=10)

_metric = (
  J=0.125,
  ξx=2.0,
  ξ̂x=16.0,
  ξy=-0.0,
  ξ̂y=-0.0,
  ηx=-0.0,
  η̂x=-0.0,
  ηy=4.0,
  η̂y=32.0,
  ξt=0.0,
  ηt=0.0,
)

cell_center_metrics_pass = true
for idx in mesh.iterators.cell.domain
  @show mesh.cell_center_metrics[idx]
  # if !(mesh.cell_center_metrics[idx] == _metric)
  #   cell_center_metrics_pass = false
  #   break
  # end
end
# @test cell_center_metrics_pass

# ilo_c, ihi_c, jlo_c, jhi_c = mesh.domain_limits.cell

# i₊½_pass = true
# for idx in mesh.iterators.cell.domain
#   @show mesh.edge_metrics.i₊½[idx]
#   # if !(mesh.edge_metrics.i₊½[idx] == _metric)
#   #   i₊½_pass = false
#   #   break
#   # end
# end

# for j in jlo_c:jhi_c
#   if !(mesh.edge_metrics.i₊½[ihi_c, j] == _metric)
#     i₊½_pass = false
#     break
#   end
# end
# @test i₊½_pass

# j₊½_pass = true
# for idx in mesh.iterators.cell.domain
#   if !(mesh.edge_metrics.j₊½[idx] == _metric)
#     j₊½_pass = false
#     break
#   end
# end

# for i in ilo_c:ihi_c
#   if !(mesh.edge_metrics.j₊½[i, jlo_c - 1] == _metric)
#     j₊½_pass = false
#     break
#   end
# end
# @test j₊½_pass

# # @test jacobian_matrix(mesh, (2, 2)) == @SMatrix [
# #   0.5 0.0
# #   0.0 0.25
# # ]

# # cell_area = 0.5 * 0.25
# # @test jacobian(mesh, (2, 3)) == cell_area

# # @test inv(jacobian_matrix(mesh, (2, 3))) == @SMatrix [
# #   2.0 0.0
# #   0.0 4.0
# # ]

# # @test inv(jacobian(mesh, (2, 3))) == 1 / cell_area

# # @benchmark conservative_metrics($mesh, (2, 3))
# # @benchmark inv(jacobian_matrix($mesh, (2, 3)))

# # bm1 = @benchmark metrics($mesh, (2, 3))
# # @test bm1.allocs == 0

# # bm2 = @benchmark jacobian_matrix($mesh, (2, 2))
# # @test bm2.allocs == 0

# # bm3 = @benchmark jacobian($mesh, (2, 2))
# # @test bm3.allocs == 0

# @test coord(mesh, (ilo_c, jlo_c)) == [0, 1]
# @test coord(mesh, (ilo_c + 1, jlo_c + 1)) == [0.5, 1.25]
# @test coord(mesh, (ilo_c + 1.5, jlo_c + 1.5)) == [0.75, 1.375]
# @test coord(mesh, (2.5, 2.5)) == centroid(mesh, (2, 2))
# @test centroid(mesh, (ilo_c, jlo_c)) == [0.25, 1.125]

# xy_coords = coords(mesh)
# centroid_coords = centroids(mesh)
# @test size(centroid_coords) == (2, 4, 8)
# @test size(xy_coords) == (2, 5, 9)
# # end

@testset "2D wavy mesh" begin
  using Test
  using Plots
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

    I₁ = (m_i₊½.ξ̂x - m_i₋½.ξ̂x) + (m_j₊½.η̂x - m_j₋½.η̂x)
    I₂ = (m_i₊½.ξ̂y - m_i₋½.ξ̂y) + (m_j₊½.η̂y - m_j₋½.η̂y)
    I₁_passes = abs(I₁) < eps()
    I₂_passes = abs(I₂) < eps()
  end
  @test I₁_passes
  @test I₂_passes
  # # cm_ij = conservative_metrics(mesh, (20.5, 20.5))
  # cm_ip1 = conservative_metrics(mesh, (21, 20))

  # cm_jm1 = conservative_metrics(mesh, (20, 20))
  # cm_jp1 = conservative_metrics(mesh, (20, 21))

  # # Jij = CurvilinearGrids.jacobian(mesh, (20, 20))
  # Jip1 = inv(CurvilinearGrids.jacobian(mesh, (21, 20)))
  # Jim1 = inv(CurvilinearGrids.jacobian(mesh, (20, 20)))
  # Jjp1 = inv(CurvilinearGrids.jacobian(mesh, (20, 21)))
  # Jjm1 = inv(CurvilinearGrids.jacobian(mesh, (20, 20)))

  # # @show Jip1, Jim1, Jjp1, Jjm1

  # m_jm1 = metrics(mesh, (20, 20))
  # m_jp1 = metrics(mesh, (20, 21))
  # m_im1 = metrics(mesh, (20, 20))
  # m_ip1 = metrics(mesh, (21, 20))
  # @show (m_ip1.ξx / Jip1), cm_ip1.ξ̂x
  # # @show m_ij.J

  # m = (m_ip1.ξx / Jip1 - m_im1.ξx / Jim1) + (m_jp1.ηx / Jjp1 - m_jm1.ηx / Jjm1)
  # cm = (cm_ip1.ξ̂x - cm_im1.ξ̂x) + (cm_jp1.η̂x - cm_jm1.η̂x)

  # @show m cm
  # (m_ip1.ξx / Jip1 - m_im1.ξx / Jim1)
  # (cm_ip1.ξ̂x - cm_im1.ξ̂x)

  # @show m_ip1.ξx, Jip1
  # @show cm_ip1.ξ̂x

  # @show cm
  # @show m
  # @show m_ip1.ξx * Jip1
  # @show m_ip1.ξx / Jip1

  # Jm = jacobian_matrix(mesh, (20, 20))
  # J⁻¹ = det(inv(Jm))
  # _J = det(Jm)
  # ξx = Jm[1, 1]
  # _ξ̂x = ξx * J⁻¹
  # _ξ̂x = ξx / _J

  # Jinv = inv(Jm)
  # # Jinv * det(Jinv)

  # cm_im1

  # CI = cell_indices(mesh)
  # total_area = 0.0
  # for idx in CI
  #   global total_area += jacobian(mesh, idx)
  # end
  # @test total_area ≈ 1.0

  # xyn = coords(mesh)

  # xn = @view xyn[1, :, :]
  # yn = @view xyn[2, :, :]

  # ∂x_∂ξ, ∂x_∂η = _fd_metrics(xn, 4)
  # ∂y_∂ξ, ∂y_∂η = _fd_metrics(yn, 4)

  # ∂x_∂ξ⁶, ∂x_∂η⁶ = _fd_metrics(xn, 6)
  # ∂y_∂ξ⁶, ∂y_∂η⁶ = _fd_metrics(yn, 6)
  # p1 = heatmap(∂x_∂ξ'; title="∂x_∂ξ")
  # p2 = heatmap(∂x_∂η'; title="∂x_∂η")
  # p3 = heatmap(∂y_∂ξ'; title="∂y_∂ξ")
  # p4 = heatmap(∂y_∂η'; title="∂y_∂η")

  # l = @layout [a b; c d]

  # pall = plot(p1, p2, p3, p4; layout=l)
  # display(pall)

  # metric_diff = ∂x_∂ξ⁶ .- ∂x_∂ξ
  # pdiff = heatmap(metric_diff; title="4th - 6th")

  # println(extrema(metric_diff))
  # pdiff
  nothing
end
