@testset "MappedGrid1D -- Rectilinear" begin
  ni = 40
  x0, x1 = (0.0, 20.0)
  Δx = (x1 - x0) / ni

  x(t, ξ, p) = x0 + (ξ - 1) * Δx
  mesh = MappedGrid(x, (;), (ni,), 5)

  cm = cell_metrics(mesh)
  fm = face_metrics(mesh)
  cell_domain = mesh.iterators.cell.domain
  @test all(cm.forward[idx].J ≈ 0.5 for idx in cell_domain)
  @test all(cm.inverse[idx][1, 1] ≈ 2.0 for idx in cell_domain)
  @test all(fm[1].conserved[idx][1, 1] ≈ 1.0 for idx in cell_domain)

  @test centroid(mesh, mesh.nhalo + 1).x ≈ 0.25
  @test collect(centroids(mesh)) == collect(0.25:0.5:19.75)
  @test collect(coords(mesh)) == collect(0:0.5:20.0)
end

@testset "DiscreteGrid1D -- Rectilinear, Halo Geometry Defined" begin
  ni = 40
  x0, x1 = (0.0, 20.0)
  x = collect(range(x0, x1, ni + 1))
  mesh = DiscreteGrid(x, 5; halo_coords_included=true)

  cm = cell_metrics(mesh)
  fm = face_metrics(mesh)
  cell_domain = mesh.iterators.cell.full
  @test all(cm.forward[idx].J ≈ 0.5 for idx in cell_domain)
  @test all(cm.inverse[idx][1, 1] ≈ 2.0 for idx in cell_domain)

  iaxis = 1
  i₊½_domain = expand(cell_domain, iaxis, -1)
  @test all(fm[1].conserved[idx][1, 1] ≈ 1.0 for idx in i₊½_domain)

  @test mesh.node_coordinates[1] == collect(0:0.5:20)
  @test mesh.centroid_coordinates[1] == collect(0.25:0.5:19.75)
  @test centroid(mesh, mesh.nhalo).x ≈ 2.25
  @test collect(centroids(mesh)) == collect(2.75:0.5:17.25)
  @test collect(coords(mesh)) == collect(2.5:0.5:17.5)
end

@testset "SphericalGrid1D" begin
  ni = 41
  x0, x1 = (0.0, 20.0)
  snap_to_axis = true
  mesh = SphericalGrid1D(range(x0, x1, ni), :meg6, snap_to_axis)

  cell_domain = mesh.iterators.cell.domain
  node_domain = mesh.iterators.node.domain
  @test all(mesh.cell_center_metrics.J[cell_domain] .≈ 0.5)
  @test all(mesh.cell_center_metrics.ξ.x₁[cell_domain] .≈ 2.0)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₁[cell_domain] .≈ 1.0)

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
  @test all(mesh.cell_center_metrics.J[cell_domain] .≈ 0.5)
  @test all(mesh.cell_center_metrics.ξ.x₁[cell_domain] .≈ 2.0)
  @test all(mesh.edge_metrics.i₊½.ξ̂.x₁[cell_domain] .≈ 1.0)

  @test centroid(mesh, mesh.nhalo + 1).x ≈ 0.25
  @test centroids(mesh) == 0.25:0.5:19.75
  @test coords(mesh) == 0:0.5:20.0

  r1 = 0.5
  @test radius(mesh, (mesh.nhalo + 2,)) ≈ r1
  @test cellvolume(mesh, (mesh.nhalo + 1,)) ≈ pi * r1^2
end
