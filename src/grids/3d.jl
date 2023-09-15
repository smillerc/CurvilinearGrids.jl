
struct CurvilinearMesh3D{T1,T2,T3,T4} <: AbstractCurvilinearMesh
  x::T1
  y::T2
  z::T3
  jacobi_func::T4
  nhalo::Int
  nnodes::NTuple{3,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi, :klo, :khi),NTuple{6,Int}}
end

function CurvilinearMesh3D(x::Function, y::Function, z::Function, (n_ξ, n_η, n_ζ), nhalo)
  jacobi_func = _setup_jacobi(x, y, z)
  nnodes = (n_ξ, n_η, n_ζ)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  nk_cells = n_ζ - 1
  lo = nhalo + 1

  limits = (
    ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo, klo=lo, khi=nk_cells - nhalo
  )
  return CurvilinearMesh3D(x, y, z, jacobi_func, nhalo, nnodes, limits)
end

function coords(m::CurvilinearMesh3D)
  xyz = zeros(3, m.nnodes...)
  for k in axes(xyz, 4)
    for j in axes(xyz, 3)
      for i in axes(xyz, 2)
        xyz[1, i, j, k] = m.x(i, j, k)
        xyz[2, i, j, k] = m.y(i, j, k)
        xyz[3, i, j, k] = m.z(i, j, k)
      end
    end
  end

  return xyz
end

function _setup_jacobi(x, y, z)

  # Make sure to define these ahead of time
  # x((ξ, η, ζ)) = x(ξ, η, ζ)
  # y((ξ, η, ζ)) = y(ξ, η, ζ)
  # z((ξ, η, ζ)) = z(ξ, η, ζ)

  ∇x(ξ, η, ζ) = gradient(x, @SVector [ξ, η, ζ]) # ∂x∂ξ, ∂x∂η, ∂x∂ζ
  ∇y(ξ, η, ζ) = gradient(y, @SVector [ξ, η, ζ]) # ∂z∂ξ, ∂z∂η, ∂z∂ζ
  ∇z(ξ, η, ζ) = gradient(z, @SVector [ξ, η, ζ]) # ∂y∂ξ, ∂y∂η, ∂y∂ζ

  function jacobi_matrix(ξ, η, ζ)
    return SMatrix{3,3,Float64}(∇x(ξ, η, ζ)..., ∇y(ξ, η, ζ)..., ∇z(ξ, η, ζ)...)
  end

  return jacobi_matrix
end

function _get_metrics(_jacobi_matrix::SMatrix)
  inv_jacobi_matrix = inv(_jacobi_matrix)
  J⁻¹ = det(inv_jacobi_matrix)

  ξx = inv_jacobi_matrix[1, 1]
  ξy = inv_jacobi_matrix[2, 1]
  ξz = inv_jacobi_matrix[3, 1]
  ηx = inv_jacobi_matrix[1, 2]
  ηy = inv_jacobi_matrix[2, 2]
  ηz = inv_jacobi_matrix[3, 2]
  ζx = inv_jacobi_matrix[1, 3]
  ζy = inv_jacobi_matrix[2, 3]
  ζz = inv_jacobi_matrix[3, 3]

  ξ̂x = ξx * J⁻¹
  ξ̂y = ξy * J⁻¹
  ξ̂z = ξz * J⁻¹
  η̂x = ηx * J⁻¹
  η̂y = ηy * J⁻¹
  η̂z = ηz * J⁻¹
  ζ̂x = ζx * J⁻¹
  ζ̂y = ζy * J⁻¹
  ζ̂z = ζz * J⁻¹

  return (ξ̂x, ξ̂y, ξ̂z, η̂x, η̂y, η̂z, ζ̂x, ζ̂y, ζ̂z)
end

metrics(m::CurvilinearMesh3D, CI::CartesianIndex) = metrics(m, CI.I...)
metrics(m::CurvilinearMesh3D, (i, j, k)) = metrics(m, i, j, k)
function metrics(m::CurvilinearMesh3D, i, j, k)
  jacobi = m.jacobi_func(i, j, k)
  return _get_metrics(jacobi)
end
