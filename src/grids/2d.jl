
"""
CurvilinearGrid2D

# Fields
 - `x`: Node function; e.g., x(i,j)
 - `y`: Node function; e.g., y(i,j)
 - `jacobian_matrix_func`: jacobian matrix, e.g., J(i,j)
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CurvilinearGrid2D{T1,T2,T3} <: AbstractCurvilinearGrid
  x::T1
  y::T2
  jacobian_matrix_func::T3
  nhalo::Int
  nnodes::NTuple{2,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi),NTuple{4,Int}}
end

function CurvilinearGrid2D(x::Function, y::Function, (n_ξ, n_η), nhalo)
  jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (n_ξ, n_η)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  lo = nhalo + 1
  limits = (ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo)

  return CurvilinearGrid2D(x, y, jacobian_matrix_func, nhalo, nnodes, limits)
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

  ∂x∂ξ(ξ, η) = ForwardDiff.derivative(ξ -> x(ξ, η), η)
  ∂x∂η(ξ, η) = ForwardDiff.derivative(η -> x(ξ, η), ξ)
  ∂y∂ξ(ξ, η) = ForwardDiff.derivative(ξ -> y(ξ, η), η)
  ∂y∂η(ξ, η) = ForwardDiff.derivative(η -> y(ξ, η), ξ)

  jacobian_matrix(ξ, η) = @SMatrix [
    ∂x∂ξ(ξ, η) ∂x∂η(ξ, η)
    ∂y∂ξ(ξ, η) ∂y∂η(ξ, η)
  ]

  # xy((i, j)) = @SVector [x((i, j)), y((i, j))]
  # jac(i, j) = ForwardDiff.jacobian(xy, @SVector [i, j])

  return jacobian_matrix
end

# Get the conservative metrics, e.g. normalized by the Jacobian
function conservative_metrics(m::CurvilinearGrid2D, (i, j))
  _jacobian_matrix = m.jacobian_matrix_func(i, j)
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J = det(_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  return (
    ξ̂x=ξx / J,
    ξ̂y=ξy / J,
    η̂x=ηx / J,
    η̂y=ηy / J,
    ξt=zero(eltype(_jacobian_matrix)),
    ηt=zero(eltype(_jacobian_matrix)),
  )
end

function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Integer})
  _jacobian_matrix = m.jacobian_matrix_func(i, j) # -> SMatrix{2,2,T}
  inv_jacobian_matrix = inv(_jacobian_matrix)

  return (
    ξx=inv_jacobian_matrix[1, 1],
    ξy=inv_jacobian_matrix[1, 2],
    ηx=inv_jacobian_matrix[2, 1],
    ηy=inv_jacobian_matrix[2, 2],
    ξt=zero(eltype(_jacobian_matrix)),
    ηt=zero(eltype(_jacobian_matrix)),
  )
end

function metrics(m::CurvilinearGrid2D, (i, j), (vx, vy))
  static = metrics(m, (i, j))
  @unpack ξx, ξy, ηx, ηy = static

  return merge(static, (ξt=-(vx * ξx + vy * ξy), ηt=-(vx * ηx + vy * ηy)))
end

"""
    coords(mesh::CurvilinearGrid2D, T=Float64) -> Array{Real}

Return the array of coordinate points, indexed as `[xy,i,j]`. 
This does _not_ include halo regions since the geometry can be undefined.
"""
function coords(m::CurvilinearGrid2D, T=Float64)
  xy = zeros(T, 2, m.nnodes...)

  @inbounds for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i, j)
      xy[2, i, j] = m.y(i, j)
    end
  end

  return xy
end

"""
    centroids(m::CurvilinearGrid2D, T=Float64) -> Array{Real}

Return the array of coordinate points, indexed as `[xy,i,j]`. 
This does _not_ include halo regions since the geometry can be undefined.
"""
function centroids(m::CurvilinearGrid2D, T=Float64)
  xy = zeros(T, 2, (m.nnodes .- 1)...)

  @inbounds for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i + 0.5, j + 0.5)
      xy[2, i, j] = m.y(i + 0.5, j + 0.5)
    end
  end

  return xy
end
