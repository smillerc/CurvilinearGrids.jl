
struct CurvilinearGrid2D{T,T1,T2,T3} <: AbstractCurvilinearGrid
  x::T1
  y::T2
  jacobian_matrix_func::T3
  nhalo::Int
  nnodes::NTuple{2,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi),NTuple{4,Int}}
  cell_center_metrics::Array{Metrics2D{T},2}
  edge_metrics::NamedTuple{
    (:ξ, :η), # one for each edge direction
    NTuple{
      2, # (ξ,η)
      Array{Metrics2D{T},2}, # an array to hold the 2d metrics
    },
  }

  J::Array{T,2} # cell-centered Jacobian J
  # J⁻¹::Array{T,2} # cell-centered inverse Jacobian J⁻¹
end

function CurvilinearGrid2D(x::Function, y::Function, (n_ξ, n_η), nhalo)
  jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (n_ξ, n_η)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  lo = nhalo + 1
  limits = (ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo)

  cell_dims = (ni_cells, nj_cells)

  cell_center_metrics = _make_empty_metric_array(cell_dims)

  # metric terms for the {i,j,k}±1/2 edges
  edge_metrics = (
    ξ=_make_empty_metric_array(cell_dims), η=_make_empty_metric_array(cell_dims)
  )

  J = zeros(cell_dims)

  return CurvilinearGrid2D(
    x, y, jacobian_matrix_func, nhalo, nnodes, limits, cell_center_metrics, edge_metrics, J
  )
end

function coords(m::CurvilinearGrid2D)
  xy = zeros(2, m.nnodes...)
  for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i, j)
      xy[2, i, j] = m.y(i, j)
    end
  end

  return xy
end

function centroids(m::CurvilinearGrid2D)
  xy = zeros(2, (m.nnodes .- 1)...)
  for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i + 0.5, j + 0.5)
      xy[2, i, j] = m.y(i + 0.5, j + 0.5)
    end
  end

  return xy
end

"""
    update_metrics(m::CurvilinearGrid3D)

Update the conservative grid metrics (ξ̂x, ξ̂y, ...). This only needs to be called once for a 
static mesh. In a dynamic mesh, this needs to be called each time the mesh moves.
"""
function update_metrics(m::CurvilinearGrid2D)
  # cell-centered metrics
  @inline for I in CartesianIndices(m.cell_center_metrics)
    ξ, η = Tuple(I)
    m.cell_center_metrics[I] = _get_metrics(m, (ξ, η))
  end

  # edge metric terms
  @inline for I in CartesianIndices(m.cell_center_metrics)
    ξ, η = Tuple(I)
    m.edge_metrics.ξ[I] = _get_metrics(m, (ξ + 0.5, η))
    m.edge_metrics.η[I] = _get_metrics(m, (ξ, η + 0.5))
  end

  return nothing
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

# Get the grid metrics for a static grid
function _get_metrics(_jacobian_matrix::SMatrix{2,2,T}) where {T}
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J = det(_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  ξ̂x = J * ξx
  ξ̂y = J * ξy
  η̂x = J * ηx
  η̂y = J * ηy

  return (ξ̂x=ξ̂x, ξ̂y=ξ̂y, η̂x=η̂x, η̂y=η̂y)
end

# Get the grid metrics for a dynamic grid -- grid velocities (vx,vy,vz) must be provided
function _get_metrics(_jacobian_matrix::SMatrix{2,2,T}, (vx, vy)) where {T}
  inv_jacobian_matrix = inv(_jacobian_matrix)
  J⁻¹ = det(inv_jacobian_matrix)

  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  ξ̂x = J⁻¹ * ξx
  ξ̂y = -J⁻¹ * ξy
  η̂x = -J⁻¹ * ηx
  η̂y = J⁻¹ * ηy

  # temporal metrics
  ξt = -(vx * ξx + vy * ξy)
  ηt = -(vx * ηx + vy * ηy)

  return (ξ̂x=ξ̂x, ξ̂y=ξ̂y, η̂x=η̂x, η̂y=η̂y, ξt=ξt, ηt=ηt)
end
