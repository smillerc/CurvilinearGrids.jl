
@testset "2D Cylindrical Mesh" begin
  function rzmesh(ni, nk)
    r0, r1 = (0, 2)
    z0, z1 = (0, 3)

    # linear distributions
    r(i, k) = r0 + (r1 - r0) * ((i - 1) / (ni - 1))
    z(i, k) = z0 + (z1 - z0) * ((k - 1) / (nk - 1))

    return (r, z)
  end

  # begin
  nr, nz = (5, 11)
  nhalo = 0
  r, z = rzmesh(nr, nz)
  snap_to_axis = true
  mesh = CylindricalGrid2D(r, z, (nr, nz), nhalo, snap_to_axis)
  # cmesh = CurvilinearGrid2D(r, z, (nr, nz), nhalo)
  # nothing
  # end

  # # J = jacobian_matrix(mesh, (2, 1))
  # jacobian_matrix(mesh, (2, 1))

  # metrics(mesh, (2, 2))
  # radius(mesh, (0, 20))

  domain = mesh.iterators.cell.domain

  I₁_passes = true
  I₂_passes = true
  for idx in domain
    i, j = idx.I .+ 0.5
    m_i₊½ = conservative_metrics(mesh, (i + 0.5, j))
    m_j₊½ = conservative_metrics(mesh, (i, j + 0.5))
    m_i₋½ = conservative_metrics(mesh, (i - 0.5, j))
    m_j₋½ = conservative_metrics(mesh, (i, j - 0.5))

    I₁ = (m_i₊½.ξ̂.x₁ - m_i₋½.ξ̂.x₁) + (m_j₊½.η̂.x₁ - m_j₋½.η̂.x₁)
    I₂ = (m_i₊½.ξ̂.x₂ - m_i₋½.ξ̂.x₂) + (m_j₊½.η̂.x₂ - m_j₋½.η̂.x₂)

    I₁_passes = abs(I₁) < eps()
    I₂_passes = abs(I₂) < eps()
    if !(I₁_passes && I₂_passes)
      break
    end
  end
  @test I₁_passes
  @test I₂_passes
end
