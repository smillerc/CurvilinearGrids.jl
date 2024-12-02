@testset "1D Mesh Rectlinear Mesh" begin
  nhalo = 3
  ni = 4
  x0, x1 = (0.0, 2.0)
  mesh = RectlinearGrid((x0, x1), ni, nhalo)

  cell_domain = mesh.iterators.cell.domain
  node_domain = mesh.iterators.node.domain
  @test all(mesh.cell_center_metrics.J[cell_domain] .≈ 0.5)
  @test all(mesh.cell_center_metrics.ξ.x₁[cell_domain] .≈ 2.0)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₁[cell_domain] .≈ 1.0)
  @test all(mesh.edge_metrics.i₊½.J[cell_domain] .≈ 0.5)

  centroid(mesh, 4).x ≈ 0.25
  @test mesh.centroid_coordinates.x[cell_domain] == [0.25, 0.75, 1.25, 1.75]
  @test mesh.node_coordinates.x[node_domain] == [0.0, 0.5, 1.0, 1.5, 2.0]
  @test centroids(mesh) == [0.25, 0.75, 1.25, 1.75]
  @test coords(mesh) == [0.0, 0.5, 1.0, 1.5, 2.0]
end
