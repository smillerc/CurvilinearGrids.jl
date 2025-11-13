@testset "2D Cylindrical Mesh" begin
  # begin
  r0, r1 = (0, 2)
  z0, z1 = (0, 4)
  nr, nz = (8, 10)
  snap_to_axis = true
  symmetry_axis = :x # rotate about the pole axis
  mesh = AxisymmetricGrid2D(
    (r0, z0), (r1, z1), (nr, nz), :meg6_symmetric, snap_to_axis, symmetry_axis
  )

  domain = mesh.iterators.cell.domain

  gcl_identities, max_vals = gcl(mesh.edge_metrics, mesh.iterators.cell.domain, eps())
  @test all(gcl_identities)

  @test centroid_radius(mesh, domain[1, 1]) == 0.2

  dx = 0.25
  dy = 0.4
  rc = 3.8
  @test centroid_radius(mesh, domain[1, end]) == rc
  @test cellvolume(mesh, domain[1, end]) ≈ dx * dy * rc * 2pi

  save_vtk(mesh, "rz_axisym")
end

@testset "Cylindrical vs Curvilinear Mesh" begin
  # begin
  r0, r1 = (0, 2)
  z0, z1 = (0, 4)
  nr, nz = (8, 10)
  snap_to_axis = true
  symmetry_axis = :x # rotate about the pole axis
  m1 = AxisymmetricGrid2D(
    (r0, z0), (r1, z1), (nr, nz), :meg6_symmetric, snap_to_axis, symmetry_axis
  )
  m2 = CurvilinearGrid2D((r0, z0), (r1, z1), (nr, nz), :meg6_symmetric)

  @test all(m1.cell_center_metrics.J .≈ m2.cell_center_metrics.J)
  @test all(m1.cell_center_metrics.ξ.x₁ .≈ m2.cell_center_metrics.ξ.x₁)
  @test all(m1.cell_center_metrics.ξ.x₂ .≈ m2.cell_center_metrics.ξ.x₂)
  @test all(m1.cell_center_metrics.ξ.t .≈ m2.cell_center_metrics.ξ.t)
  @test all(m1.cell_center_metrics.η.x₁ .≈ m2.cell_center_metrics.η.x₁)
  @test all(m1.cell_center_metrics.η.x₂ .≈ m2.cell_center_metrics.η.x₂)
  @test all(m1.cell_center_metrics.η.t .≈ m2.cell_center_metrics.η.t)
  @test all(m1.cell_center_metrics.ξ̂.x₁ .≈ m2.cell_center_metrics.ξ̂.x₁)
  @test all(m1.cell_center_metrics.ξ̂.x₂ .≈ m2.cell_center_metrics.ξ̂.x₂)
  @test all(m1.cell_center_metrics.ξ̂.t .≈ m2.cell_center_metrics.ξ̂.t)
  @test all(m1.cell_center_metrics.η̂.x₁ .≈ m2.cell_center_metrics.η̂.x₁)
  @test all(m1.cell_center_metrics.η̂.x₂ .≈ m2.cell_center_metrics.η̂.x₂)
  @test all(m1.cell_center_metrics.η̂.t .≈ m2.cell_center_metrics.η̂.t)
  @test all(m1.cell_center_metrics.x₁.ξ .≈ m2.cell_center_metrics.x₁.ξ)
  @test all(m1.cell_center_metrics.x₁.η .≈ m2.cell_center_metrics.x₁.η)
  @test all(m1.cell_center_metrics.x₂.ξ .≈ m2.cell_center_metrics.x₂.ξ)
  @test all(m1.cell_center_metrics.x₂.η .≈ m2.cell_center_metrics.x₂.η)
end
