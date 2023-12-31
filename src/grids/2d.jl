
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
  dim = 2
  check_nargs(x, dim, :x)
  check_nargs(y, dim, :y)
  test_coord_func(x, dim, :x)
  test_coord_func(y, dim, :y)

  coord(i, j) = @SVector [x(i, j), y(i, j)]
  jacobian_matrix_func(i, j) = ForwardDiff.jacobian(x -> coord(x[1], x[2]), @SVector [i, j])

  # jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (n_ξ, n_η)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  lo = nhalo + 1
  limits = (ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo)

  return CurvilinearGrid2D(x, y, jacobian_matrix_func, nhalo, nnodes, limits)
end

# Get the conservative metrics, e.g. normalized by the Jacobian
@inline function conservative_metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
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

@inline function conservative_metrics(
  m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, (vx, vy)
)
  static = conservative_metrics(m, (i, j))
  @unpack ξ̂x, ξ̂y, η̂x, η̂y = static

  return merge(static, (
    ξt=-(vx * ξ̂x + vy * ξ̂y), # dynamic / moving mesh terms
    ηt=-(vx * η̂x + vy * η̂y), # dynamic / moving mesh terms
  ))
end

@inline function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo) # -> SMatrix{2,2,T}
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

@inline function cell_metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return metrics(m, (i + 0.5, j + 0.5))
end

@inline function metrics_with_jacobian(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  _jacobian_matrix = checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo)) # -> SMatrix{2,2,T}
  inv_jacobian_matrix = inv(_jacobian_matrix)

  return (
    ξx=inv_jacobian_matrix[1, 1],
    ξy=inv_jacobian_matrix[1, 2],
    ηx=inv_jacobian_matrix[2, 1],
    ηy=inv_jacobian_matrix[2, 2],
    ξt=zero(eltype(_jacobian_matrix)),
    ηt=zero(eltype(_jacobian_matrix)),
    J=det(_jacobian_matrix),
  )
end

@inline function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, (vx, vy))
  static = metrics(m, (i, j))
  @unpack ξx, ξy, ηx, ηy = static

  return merge(static, (
    ξt=-(vx * ξx + vy * ξy), # dynamic / moving mesh terms
    ηt=-(vx * ηx + vy * ηy), # dynamic / moving mesh terms
  ))
end

function jacobian_matrix(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return checkeps(m.jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
end

function jacobian(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return det(jacobian_matrix(m, (i, j)))
end

"""
    coords(mesh::CurvilinearGrid2D, T=Float64) -> Array{Real}

Return the array of coordinate points, indexed as `[xy,i,j]`.
This does _not_ include halo regions since the geometry can be undefined.
"""
function coords!(xy::Array{T,3}, m::CurvilinearGrid2D) where {T}
  @inbounds for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i, j)
      xy[2, i, j] = m.y(i, j)
    end
  end

  return nothing
end

function coords(m::CurvilinearGrid2D, T=Float64)
  xy = zeros(T, 2, m.nnodes...)
  coords!(xy, m)
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
