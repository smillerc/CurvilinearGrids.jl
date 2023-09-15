
struct CurvilinearMesh2D{T1,T2,T3} <: AbstractCurvilinearMesh
  x::T1
  y::T2
  jacobian_matrix_func::T3
  nhalo::Int
  nnodes::NTuple{2,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi),NTuple{4,Int}}
end

function CurvilinearMesh2D(x::Function, y::Function, (n_ξ, n_η), nhalo)
  jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (n_ξ, n_η)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  lo = nhalo + 1
  limits = (ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo)
  return CurvilinearMesh2D(x, y, jacobian_matrix_func, nhalo, nnodes, limits)
end

function coords(m::CurvilinearMesh2D)
  xy = zeros(2, m.nnodes...)
  for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i, j)
      xy[2, i, j] = m.y(i, j)
    end
  end

  return xy
end

function centroids(m::CurvilinearMesh2D)
  xy = zeros(2, (m.nnodes .- 1)...)
  for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i + 0.5, j + 0.5)
      xy[2, i, j] = m.y(i + 0.5, j + 0.5)
    end
  end

  return xy
end

function _setup_jacobian_func(x, y)

  # for the gradient calls to work, the x,y,z functions
  # need to work with a single argument, so we make these
  # functions below:
  error_found = false
  for func in (x, y)
    try
      func((1, 1))
    catch
      @warn(
        "The $(Symbol(func)) function needs to have definitions for both $(Symbol(func))(i,j) and $(Symbol(func))((i,j)), or else the gradient calculations will fail",
      )
      error_found = true
    end
  end

  if error_found
    error("Grid function definitions aren't complete... see the warning messages")
  end

  ∂x∂ξ(ξ, η) = derivative(ξ -> x(ξ, η), η)
  ∂x∂η(ξ, η) = derivative(η -> x(ξ, η), ξ)
  ∂y∂ξ(ξ, η) = derivative(ξ -> y(ξ, η), η)
  ∂y∂η(ξ, η) = derivative(η -> y(ξ, η), ξ)

  jacobian_matrix(ξ, η) = @SMatrix [
    ∂x∂ξ(ξ, η) ∂x∂η(ξ, η)
    ∂y∂ξ(ξ, η) ∂y∂η(ξ, η)
  ]

  return jacobian_matrix
end

function _get_metrics(_jacobian_matrix::SMatrix{2,2,T}) where {T}
  J = det(_jacobian_matrix)

  ∂x∂ξ = _jacobian_matrix[1, 1]
  ∂x∂η = _jacobian_matrix[1, 2]
  ∂y∂ξ = _jacobian_matrix[2, 1]
  ∂y∂η = _jacobian_matrix[2, 2]

  return (ξ̂x=J * ∂y∂η, ξ̂y=-J * ∂x∂η, η̂x=-J * ∂y∂ξ, η̂y=J * ∂x∂ξ)
end
