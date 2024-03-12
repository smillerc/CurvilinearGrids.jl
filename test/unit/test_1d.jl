
@testset "1D Mesh Rectlinear Mesh" begin
  function getmesh()
    nhalo = 3
    ni = 5
    x0, x1 = (0.0, 2.0)
    x(i) = (x0 + (x1 - x0) * ((i - 1) / (ni - 1)))
    return CurvilinearGrid1D(x, ni, nhalo)
  end

  m = getmesh()

  m.cell_center_metrics
  m.edge_metrics

  @test all(m.cell_center_metrics.J .≈ 0.5)
  @test all(m.cell_center_metrics.ξ.x₁ .≈ 2.0)
  @test all(m.edge_metrics.i₊½.ξ̂.x₁ .≈ 1.0)
  @test all(m.edge_metrics.i₊½.J .≈ 0.5)

  centroid(m, 4) ≈ 0.25
  @test m.centroid_coordinates.x ==
    [-1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25]

  @test m.node_coordinates.x == [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

  # @test metrics(m, 1) == (ξx=2.0, ξt=0.0)
  # @test conservative_metrics(m, 1) == (ξ̂x=4.0, ξt=0.0)
  # @test jacobian_matrix(m, 1) == @SMatrix [0.5]
  # @test jacobian(m, 1) == 0.5
  # @test inv(jacobian_matrix(m, 1)) == @SMatrix [2.0]
  # @test inv(jacobian(m, 1)) == 2.0

  # @test metrics(m, (1,)) == (ξx=2.0, ξt=0.0)
  # @test conservative_metrics(m, (1,)) == (ξ̂x=4.0, ξt=0.0)
  # @test jacobian_matrix(m, (1,)) == @SMatrix [0.5]
  # @test jacobian(m, (1,)) == 0.5
  # @test inv(jacobian_matrix(m, (1,))) == @SMatrix [2.0]
  # @test inv(jacobian(m, (1,))) == 2.0

  # @test centroids(m) == [0.25, 0.75, 1.25, 1.75]
  # @test coords(m) == [0.0, 0.5, 1.0, 1.5, 2.0]

  bm1 = @benchmark metrics($m, 1, 0)
  @test bm1.allocs == 0

  bm2 = @benchmark conservative_metrics($m, 1, 0)
  @test bm2.allocs == 0

  # bm3 = @benchmark jacobian_matrix($m, 1, 0)
  # @test bm3.allocs == 0

  # bm4 = @benchmark jacobian($m, 1, 0)
  # @test bm4.allocs == 0

  # bm5 = @benchmark inv(jacobian_matrix($m, 1))
  # @test bm5.allocs == 0

  # bm6 = @benchmark inv(jacobian($m, 1))
  # @test bm6.allocs == 0
end
