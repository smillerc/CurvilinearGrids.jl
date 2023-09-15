
struct CurvilinearMesh3D{T1,T2,T3,T4} <: AbstractCurvilinearMesh
  x::T1
  y::T2
  z::T3
  jacobian_matrix_func::T4
  nhalo::Int
  nnodes::NTuple{3,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi, :klo, :khi),NTuple{6,Int}}
end

function CurvilinearMesh3D(x::Function, y::Function, z::Function, (n_ξ, n_η, n_ζ), nhalo)
  jacobian_matrix_func = _setup_jacobian_func(x, y, z)
  nnodes = (n_ξ, n_η, n_ζ)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  nk_cells = n_ζ - 1
  lo = nhalo + 1

  limits = (
    ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo, klo=lo, khi=nk_cells - nhalo
  )
  return CurvilinearMesh3D(x, y, z, jacobian_matrix_func, nhalo, nnodes, limits)
end

function coords(m::CurvilinearMesh3D)
  dims = m.nnodes
  xyz = zeros(3, dims...)
  @inbounds for I in CartesianIndices(dims)
    i, j, k = Tuple(I)
    xyz[1, i, j, k] = m.x(i, j, k)
    xyz[2, i, j, k] = m.y(i, j, k)
    xyz[3, i, j, k] = m.z(i, j, k)
  end

  return xyz
end

function centroids(m::CurvilinearMesh3D)
  dims = (m.nnodes .- 1)
  xyz = zeros(3, dims...)
  @inbounds for I in CartesianIndices(dims)
    i, j, k = Tuple(I)
    xyz[1, i, j, k] = m.x(i + 0.5, j + 0.5, k + 0.5)
    xyz[2, i, j, k] = m.y(i + 0.5, j + 0.5, k + 0.5)
    xyz[3, i, j, k] = m.z(i + 0.5, j + 0.5, k + 0.5)
  end

  return xyz
end

function _setup_jacobian_func(x, y, z)

  # for the gradient calls to work, the x,y,z functions
  # need to work with a single argument, so we make these
  # functions below:
  error_found = false
  for func in (x, y, z)
    try
      func((1, 1, 1))
    catch
      @warn(
        "The $(Symbol(func)) function needs to have definitions for both $(Symbol(func))(i,j,k) and $(Symbol(func))((i,j,k)), or else the gradient calculations will fail",
      )
      error_found = true
    end
  end

  if error_found
    error("Grid function definitions aren't complete... see the warning messages")
  end

  ∇x(ξ, η, ζ) = gradient(x, @SVector [ξ, η, ζ]) # ∂x∂ξ, ∂x∂η, ∂x∂ζ
  ∇y(ξ, η, ζ) = gradient(y, @SVector [ξ, η, ζ]) # ∂z∂ξ, ∂z∂η, ∂z∂ζ
  ∇z(ξ, η, ζ) = gradient(z, @SVector [ξ, η, ζ]) # ∂y∂ξ, ∂y∂η, ∂y∂ζ

  function jacobian_matrix(ξ, η, ζ)
    return SMatrix{3,3,Float64}(∇x(ξ, η, ζ)..., ∇y(ξ, η, ζ)..., ∇z(ξ, η, ζ)...)
  end

  return jacobian_matrix
end

# Get the grid metrics for a static grid
function _get_metrics(_jacobian_matrix::SMatrix{3,3,T}) where {T}
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J⁻¹ = det(inv_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ξz = inv_jacobian_matrix[1, 3]

  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]
  ηz = inv_jacobian_matrix[2, 3]

  ζx = inv_jacobian_matrix[3, 1]
  ζy = inv_jacobian_matrix[3, 2]
  ζz = inv_jacobian_matrix[3, 3]

  return (
    ξ̂x=ξx * J⁻¹,
    ξ̂y=ξy * J⁻¹,
    ξ̂z=ξz * J⁻¹,
    η̂x=ηx * J⁻¹,
    η̂y=ηy * J⁻¹,
    η̂z=ηz * J⁻¹,
    ζ̂x=ζx * J⁻¹,
    ζ̂y=ζy * J⁻¹,
    ζ̂z=ζz * J⁻¹,
  )
end

# Get the grid metrics for a dynamic grid -- grid velocities (vx,vy,vz) must be provided
function _get_metrics(_jacobian_matrix::SMatrix{3,3,T}, (vx, vy, vz)) where {T}
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J⁻¹ = det(inv_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ξz = inv_jacobian_matrix[1, 3]

  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]
  ηz = inv_jacobian_matrix[2, 3]

  ζx = inv_jacobian_matrix[3, 1]
  ζy = inv_jacobian_matrix[3, 2]
  ζz = inv_jacobian_matrix[3, 3]

  # temporal metrics
  ξt = -(vx * ξx + vy * ξy + vz * ξz)
  ηt = -(vx * ηx + vy * ηy + vz * ηz)
  ζt = -(vx * ζx + vy * ζy + vz * ζz)

  ξ̂x = ξx * J⁻¹
  ξ̂y = ξy * J⁻¹
  ξ̂z = ξz * J⁻¹
  η̂x = ηx * J⁻¹
  η̂y = ηy * J⁻¹
  η̂z = ηz * J⁻¹
  ζ̂x = ζx * J⁻¹
  ζ̂y = ζy * J⁻¹
  ζ̂z = ζz * J⁻¹

  return (
    ξ̂x=ξ̂x,
    ξ̂y=ξ̂y,
    ξ̂z=ξ̂z,
    η̂x=η̂x,
    η̂y=η̂y,
    η̂z=η̂z,
    ζ̂x=ζ̂x,
    ζ̂y=ζ̂y,
    ζ̂z=ζ̂z,
    ξt=ξt,
    ηt=ηt,
    ζt=ζt,
  )
end
