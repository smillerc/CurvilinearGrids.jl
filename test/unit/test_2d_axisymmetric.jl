@testset "2D Cylindrical Mesh" begin
  # begin
  r0, r1 = (0, 2)
  z0, z1 = (0, 4)
  nr, nz = (8, 10)
  snap_to_axis = true
  symmetry_axis = :x # rotate about the pole axis
  mesh = axisymmetric_rectilinear_grid(
    (r0, z0), (r1, z1), (nr, nz), :meg6_symmetric, snap_to_axis, symmetry_axis
  )

  domain = mesh.iterators.cell.domain

  # function gcl(mesh)
  #   I₁_passes = true
  #   I₂_passes = true

  #   ϵ = 5e-14
  #   em = mesh.edge_metrics.inverse_normalized

  #   for idx in mesh.iterators.cell.domain
  #     i, j = idx.I
  #     I₁ = (
  #       (em.i₊½.ξ̂.x₁[i, j] - em.i₊½.ξ̂.x₁[i - 1, j]) +
  #       (em.j₊½.η̂.x₁[i, j] - em.j₊½.η̂.x₁[i, j - 1])
  #     )
  #     I₂ = (
  #       (em.i₊½.ξ̂.x₂[i, j] - em.i₊½.ξ̂.x₂[i - 1, j]) +
  #       (em.j₊½.η̂.x₂[i, j] - em.j₊½.η̂.x₂[i, j - 1])
  #     )

  #     I₁_passes = abs(I₁) < ϵ
  #     I₂_passes = abs(I₂) < ϵ
  #     if !(I₁_passes && I₂_passes)
  #       @show I₁ I₂
  #       break
  #     end
  #   end
  #   @test I₁_passes
  #   @test I₂_passes
  # end

  gcl(mesh)
  @test centroid_radius(mesh, domain[1, 1]) == 0.2

  dx = 0.25
  dy = 0.4
  rc = 3.8
  @test centroid_radius(mesh, domain[1, end]) == rc
  @test cellvolume(mesh, domain[1, end]) ≈ dx * dy * rc * 2pi

  save_vtk(mesh, "rz_axisym")
end
