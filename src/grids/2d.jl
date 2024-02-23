
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
struct CurvilinearGrid2D{CO,CE,NV,EM,CM,DL,CI,CF,JF} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{2,Int}
  domain_limits::DL
  iterators::CI
  _coordinate_funcs::CF
  _jacobian_matrix_func::JF
end

function CurvilinearGrid2D(
  x::Function, y::Function, (n_ξ, n_η), nhalo; T=Float64, backend=CPU()
)
  dim = 2
  check_nargs(x, dim, :x)
  check_nargs(y, dim, :y)
  test_coord_func(x, dim, :x)
  test_coord_func(y, dim, :y)

  xy(i, j) = @SVector [x(i, j), y(i, j)]
  jacobian_matrix_func(i, j) = ForwardDiff.jacobian(x -> xy(x[1], x[2]), @SVector [i, j])

  # jacobian_matrix_func = _setup_jacobian_func(x, y)
  nnodes = (n_ξ, n_η)
  ncells = nnodes .- 1
  ni_cells, nj_cells = ncells
  lo = nhalo + 1
  limits = (
    node=(ilo=lo, ihi=n_ξ + nhalo, jlo=lo, jhi=n_η + nhalo),
    cell=(ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo),
  )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)
  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x=KernelAbstractions.zeros(backend, T, celldims),
      y=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    η=StructArray((
      x=KernelAbstractions.zeros(backend, T, celldims),
      y=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
    j₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x=KernelAbstractions.zeros(backend, T, celldims),
        y=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
  )

  coordinate_funcs = (; x, y)
  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
  ))
  _centroid_coordinates!(centroids, coordinate_funcs, domain_iterators.cell.full, nhalo)

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))
  _node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  m = CurvilinearGrid2D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    domain_iterators,
    coordinate_funcs,
    jacobian_matrix_func,
  )

  update_metrics!(m)
  # check_for_invalid_metrics(m)
  return m
end

function update_metrics!(m::CurvilinearGrid2D, t=0)

  # cell metrics
  @inbounds for idx in m.iterators.cell.full
    cell_idx = idx.I .+ 0.5
    # @unpack J, ξ, η, x, y = metrics(m, cell_idx, t)
    @unpack J, ξ, η = metrics(m, cell_idx, t)

    m.cell_center_metrics.ξ.x[idx] = ξ.x
    m.cell_center_metrics.ξ.y[idx] = ξ.y
    m.cell_center_metrics.ξ.t[idx] = ξ.t
    m.cell_center_metrics.η.x[idx] = η.x
    m.cell_center_metrics.η.y[idx] = η.y
    m.cell_center_metrics.η.t[idx] = η.t

    # m.cell_center_inv_metrics.xξ[idx] = x.ξ
    # m.cell_center_inv_metrics.yξ[idx] = y.ξ
    # m.cell_center_inv_metrics.xη[idx] = x.η
    # m.cell_center_inv_metrics.yη[idx] = y.η

    m.cell_center_metrics.J[idx] = J
  end

  # i₊½ conserved metrics
  @inbounds for idx in m.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i₊½, j)
    @unpack ξ̂, η̂, J = conservative_metrics(m, (i + 1 / 2, j), t)

    m.edge_metrics.i₊½.ξ̂.x[idx] = ξ̂.x
    m.edge_metrics.i₊½.ξ̂.y[idx] = ξ̂.y
    m.edge_metrics.i₊½.ξ̂.t[idx] = ξ̂.t
    m.edge_metrics.i₊½.η̂.x[idx] = η̂.x
    m.edge_metrics.i₊½.η̂.y[idx] = η̂.y
    m.edge_metrics.i₊½.η̂.t[idx] = η̂.t
    m.edge_metrics.i₊½.J[idx] = J
  end

  # j₊½ conserved metrics
  @inbounds for idx in m.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i, j₊½)
    @unpack ξ̂, η̂, J = conservative_metrics(m, (i, j + 1 / 2), t)

    m.edge_metrics.j₊½.ξ̂.x[idx] = ξ̂.x
    m.edge_metrics.j₊½.ξ̂.y[idx] = ξ̂.y
    m.edge_metrics.j₊½.ξ̂.t[idx] = ξ̂.t
    m.edge_metrics.j₊½.η̂.x[idx] = η̂.x
    m.edge_metrics.j₊½.η̂.y[idx] = η̂.y
    m.edge_metrics.j₊½.η̂.t[idx] = η̂.t
    m.edge_metrics.j₊½.J[idx] = J
  end

  return nothing
end

# ------------------------------------------------------------------
# Grid Metrics
# ------------------------------------------------------------------

@inline function metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real)
  _jacobian_matrix = checkeps(m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  J = det(_jacobian_matrix)

  vx, vy = grid_velocities(m, (i, j), t)
  ξt = -(vx * ξx + vy * ξy)
  ηt = -(vx * ηx + vy * ηy)

  ξ = Metric2D(ξx, ξy, ξt)
  η = Metric2D(ηx, ηy, ηt)

  return (; ξ, η, J)
end

# ------------------------------------------------------------------
# Conservative Grid Metrics; e.g. ξ̂x = ξx * J
# ------------------------------------------------------------------

@inline conservative_metrics(m::CurvilinearGrid2D, idx) = conservative_metrics(m, idx, 0)

@inline function conservative_metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real)
  _jacobian_matrix = checkeps(m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  J = det(_jacobian_matrix)

  vx, vy = grid_velocities(m, (i, j), t)
  ξt = -(vx * ξx + vy * ξy)
  ηt = -(vx * ηx + vy * ηy)

  ξ̂ = Metric2D(ξx * J, ξy * J, ξt * J)
  η̂ = Metric2D(ηx * J, ηy * J, ηt * J)

  return (; ξ̂, η̂, J)
end

# ------------------------------------------------------------------
# Jacobian related functions
# ------------------------------------------------------------------
function jacobian_matrix(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return checkeps(m._jacobian_matrix_func(i - m.nhalo, j - m.nhalo))
end

function jacobian(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return det(jacobian_matrix(m, (i, j)))
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline centroid_velocities(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline node_velocities(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

function _node_coordinates!(
  coordinates::StructArray{T,2}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the node coordinates
  @inbounds for idx in domain
    cell_idx = @. idx.I - nhalo
    coordinates.x[idx] = coordinate_functions.x(cell_idx...)
    coordinates.y[idx] = coordinate_functions.y(cell_idx...)
  end

  return nothing
end

function _centroid_coordinates!(
  centroids::StructArray{T,2}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the centroid coordinates
  @inbounds for idx in domain
    cell_idx = @. idx.I - nhalo + 0.5
    centroids.x[idx] = coordinate_functions.x(cell_idx...)
    centroids.y[idx] = coordinate_functions.y(cell_idx...)
  end

  return nothing
end
