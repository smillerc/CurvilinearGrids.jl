
using CurvilinearGrids

function rzmesh(ni, nk)
  r0, r1 = (0, 2)
  z0, z1 = (0, 3)

  # linear distributions
  r(i, k) = r0 + (r1 - r0) * ((i - 1) / (ni - 1))
  z(i, k) = z0 + (z1 - z0) * ((k - 1) / (nk - 1))

  return (r, z)
end

begin
  nr, nz = (5, 11)
  nhalo = 0
  r, z = rzmesh(nr, nz)
  mesh = CylindricalGrid2D(r, z, (nr, nz), nhalo)
  # cmesh = CurvilinearGrid2D(r, z, (nr, nz), nhalo)
  nothing
end

# J = jacobian_matrix(mesh, (2, 1))
jacobian_matrix(mesh, (2, 1))

metrics(mesh, (2, 2))
radius(mesh, (0, 20))

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

# inv(jacobian_matrix(mesh, (1, 2)))
# inv(jacobian_matrix(mesh, (2, 1, 2)))

# jacobian(mesh, (1, 2))
# jacobian(mesh, (1, 1, 2))

# metrics(mesh, (2, 2, 2), 0)

# # Extract the 2D metrics, for use in a 2d axisymmetric fluid simulation. These are as if
# # the mesh is cartesian. The axisymmetry is enforced via a source term later...
# o = 1 + eps()

# metrics(mesh, (1, 2), 0)
# metrics(mesh, (o, 2), 0)
# metrics(mesh, (2, 2), 0)
# metrics(mesh, (3, 2), 0)

# metrics(cmesh, (1, 2), 0)
# metrics(cmesh, (2, 2), 0)
# metrics(cmesh, (3, 2), 0)

# metrics(mesh, (1, 2))
# metrics(mesh, (1 + eps(), 1, 2))
# metrics(mesh, (1, 1, 2))
# metrics(mesh, (1, 1, 2))

# # How do we distinguish between 2d metrics in a planar sense or 2d in an axisymmetric case???

# # For the sake of this package's API, it needs to be _consistent_ with the geometry!

# function covariant_basis_vectors(m, (i, j))
#   J = jacobian_matrix(m, (i, j))
#   e₁₂ = @SVector [J[:, 1], J[:, 2]]
#   return e₁₂
# end

# function covariant_basis_vectors(m, (i, j, k))
#   J = jacobian_matrix(m, (i, j, k))
#   e₁₂₃ = @SVector [J[:, 1], J[:, 2], J[:, 3]]
#   return e₁₂₃
# end

# function contravariant_basis_vectors(m, (i, j))
#   J⁻¹ = inv(jacobian_matrix(m, (i, j)))
#   e¹² = @SVector [J⁻¹[1, :], J⁻¹[2, :]]
#   return e¹²
# end

# covariant_basis_vectors(mesh, (2, 2))
# covariant_basis_vectors(mesh, (2, 0.25, 3))
# contravariant_basis_vectors(mesh, (2, 2))

# coord(mesh, (1, 2, 3))
# coord(mesh, (1, 3))
