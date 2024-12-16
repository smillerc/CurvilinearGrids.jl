@testset "CylindricalGrid1D -- Rectlinear" begin
  ni = 40
  x0, x1 = (0.0, 20.0)
  mesh = rectlinear_grid((x0, x1), ni, :meg6)

  cell_domain = mesh.iterators.cell.domain
  node_domain = mesh.iterators.node.domain
  @test all(mesh.cell_center_metrics.forward.J[cell_domain] .≈ 0.5)
  @test all(mesh.cell_center_metrics.inverse.ξ.x₁[cell_domain] .≈ 2.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₁[cell_domain] .≈ 1.0)

  @test centroid(mesh, mesh.nhalo + 1).x ≈ 0.25
  @test centroids(mesh) == 0.25:0.5:19.75
  @test coords(mesh) == 0:0.5:20.0
end

@testset "SphericalGrid1D" begin
  ni = 41
  x0, x1 = (0.0, 20.0)
  snap_to_axis = true
  mesh = SphericalGrid1D(range(x0, x1, ni), :meg6, snap_to_axis)

  cell_domain = mesh.iterators.cell.domain
  node_domain = mesh.iterators.node.domain
  @test all(mesh.cell_center_metrics.forward.J[cell_domain] .≈ 0.5)
  @test all(mesh.cell_center_metrics.inverse.ξ.x₁[cell_domain] .≈ 2.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₁[cell_domain] .≈ 1.0)

  @test centroid(mesh, mesh.nhalo + 1).x ≈ 0.25
  @test centroids(mesh) == 0.25:0.5:19.75
  @test coords(mesh) == 0:0.5:20.0

  r1 = 0.5
  @test radius(mesh, (mesh.nhalo + 2,)) ≈ r1
  @test cellvolume(mesh, (mesh.nhalo + 1,)) ≈ (4 / 3) * pi * r1^3
end

@testset "CylindricalGrid1D" begin
  ni = 41
  x0, x1 = (0.0, 20.0)
  snap_to_axis = true
  mesh = CylindricalGrid1D(range(x0, x1, ni), :meg6, snap_to_axis)

  cell_domain = mesh.iterators.cell.domain
  node_domain = mesh.iterators.node.domain
  @test all(mesh.cell_center_metrics.forward.J[cell_domain] .≈ 0.5)
  @test all(mesh.cell_center_metrics.inverse.ξ.x₁[cell_domain] .≈ 2.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₁[cell_domain] .≈ 1.0)

  @test centroid(mesh, mesh.nhalo + 1).x ≈ 0.25
  @test centroids(mesh) == 0.25:0.5:19.75
  @test coords(mesh) == 0:0.5:20.0

  r1 = 0.5
  @test radius(mesh, (mesh.nhalo + 2,)) ≈ r1
  @test cellvolume(mesh, (mesh.nhalo + 1,)) ≈ pi * r1^2
end
