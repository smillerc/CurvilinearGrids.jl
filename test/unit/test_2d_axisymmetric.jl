@testset "2D Cylindrical Mesh" begin

  # begin
  r0, r1 = (0, 2)
  z0, z1 = (0, 3)
  nr, nz = (4, 10)
  nhalo = 4
  snap_to_axis = true
  symmetry_axis = :x # rotate about the pole axis
  mesh = AxisymmetricRectlinearGrid(
    (r0, z0), (r1, z1), (nr, nz), nhalo, snap_to_axis, symmetry_axis
  )

  domain = mesh.iterators.cell.domain

  function gcl(mesh, domain)
    ϵ = 5e-14
    I₁_passes = true
    I₂_passes = true
    for idx in domain
      i, j = idx.I

      ξ̂_i₊½ = mesh.edge_metrics.i₊½.ξ̂[i, j]
      ξ̂_i₋½ = mesh.edge_metrics.i₊½.ξ̂[i - 1, j]
      η̂_j₊½ = mesh.edge_metrics.j₊½.η̂[i, j]
      η̂_j₋½ = mesh.edge_metrics.j₊½.η̂[i, j - 1]

      I₁ = (ξ̂_i₊½.x₁ - ξ̂_i₋½.x₁) + (η̂_j₊½.x₁ - η̂_j₋½.x₁)
      I₂ = (ξ̂_i₊½.x₂ - ξ̂_i₋½.x₂) + (η̂_j₊½.x₂ - η̂_j₋½.x₂)

      I₁ = I₁ * (abs(I₁) >= ϵ)
      I₂ = I₂ * (abs(I₂) >= ϵ)

      # @show I₁, I₂
      I₁_passes = abs(I₁) < eps()
      I₂_passes = abs(I₂) < eps()
      if !(I₁_passes && I₂_passes)
        break
      end
    end
    @test I₁_passes
    @test I₂_passes
  end

  gcl(mesh, mesh.iterators.cell.domain)
  @test centroid_radius(mesh, domain[1, 1]) == 0.15

  dx = 0.5
  dy = 0.3
  rc = 2.85
  @test centroid_radius(mesh, domain[1, end]) == rc
  @test cellvolume(mesh, domain[1, end]) ≈ dx * dy * rc * 2pi

  save_vtk(mesh, "rz_axisym")
end
