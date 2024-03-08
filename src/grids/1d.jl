
"""
CurvilinearGrid1D

# Fields
 - `x`: Node function, e.g. x -> f(ξ)
 - `∂x∂ξ`: Derivative of x wrt ξ; ∂x∂ξ(ξ)
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CurvilinearGrid1D{CO,CE,NV,EM,CM,DL,CI,CF,DX} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  _coordinate_funcs::CF
  _∂x∂ξ::DX
end

"""
    CurvilinearGrid1D(x::Function, (n_ξ,), nhalo)

Create a `CurvilinearGrid1D` with a function `x(ξ)` and `nhalo` halo cells. `n_ξ` is the
total number of nodes/vertices (not including halo).
"""
# function CurvilinearGrid1D(x::F1, (n_ξ,), nhalo) where {F1<:Function}
#   ∂x∂ξ(ξ) = ForwardDiff.derivative(x, ξ)
#   F2 = typeof(∂x∂ξ)
#   nnodes = (n_ξ,)
#   ni_cells = n_ξ - 1

#   # limits used for looping to avoid using the halo regions
#   limits = (
#     ilo=nhalo + 1, # starting index
#     ihi=ni_cells + nhalo, # ending index
#   )

#   return CurvilinearGrid1D{F1,F2}(x, ∂x∂ξ, nhalo, nnodes, limits)
# end

function CurvilinearGrid1D(x::Function, ni::Int, nhalo; T=Float64, backend=CPU())
  dim = 1
  check_nargs(x, dim, :x)
  test_coord_func(x, dim, :x)

  function ∂x∂ξ(ξ, t)
    return ForwardDiff.derivative(x, ξ)
  end

  nnodes = ni
  ncells = nnodes - 1
  ilo = nhalo + 1
  limits = (node=(ilo=ilo, ihi=ni + nhalo), cell=(ilo=ilo, ihi=ncells + nhalo))

  nodeCI = CartesianIndices((nnodes + 2nhalo,))
  cellCI = CartesianIndices((ncells + 2nhalo,))

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coordinate_funcs = (; x)
  centroids = StructArray((x=KernelAbstractions.zeros(backend, T, celldims),))
  _centroid_coordinates!(centroids, coordinate_funcs, domain_iterators.cell.full, nhalo)

  coords = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))
  _node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  m = CurvilinearGrid1D(
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
    ∂x∂ξ,
  )

  update_metrics!(m)
  # check_for_invalid_metrics(m)
  return m
end

function update_metrics!(m::CurvilinearGrid1D, t::Real=0)
  # cell metrics
  @inbounds for idx in m.iterators.cell.full
    cell_idx, = idx.I .+ 0.5
    # @unpack J, ξ, x = metrics(m, cell_idx, t)
    @unpack J, ξ = metrics(m, cell_idx, t)
    m.cell_center_metrics.ξ.x₁[idx] = ξ.x₁
    # m.cell_center_inv_metrics.xξ[idx] = x.ξ
    m.cell_center_metrics.J[idx] = J
  end

  # get the conserved metrics at (i₊½)
  @inbounds for idx in m.iterators.cell.full
    i, = idx.I .+ 0.5 # centroid index

    @unpack ξ̂, J = conservative_metrics(m, i + 1 / 2, t)

    m.edge_metrics.i₊½.ξ̂.x₁[idx] = ξ̂.x₁
    m.edge_metrics.i₊½.J[idx] = J
  end

  return nothing
end

# ------------------------------------------------------------------
# Grid Metrics
# ------------------------------------------------------------------

@inline function metrics(m::CurvilinearGrid1D, i::Real, t::Real=0)
  _jacobian_matrix = checkeps(m._∂x∂ξ(i - m.nhalo, t))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1]
  J = det(_jacobian_matrix)

  vx = grid_velocities(m, i, t)
  ξt = -(vx * ξx)

  ξ = Metric1D(ξx, ξt)

  return (; ξ, J)
end

# ------------------------------------------------------------------
# Conservative Grid Metrics; e.g. ξ̂x = ξx * J
# ------------------------------------------------------------------

@inline function conservative_metrics(m::CurvilinearGrid1D, i::Real, t::Real=0)
  _jacobian_matrix = checkeps(m._∂x∂ξ(i - m.nhalo, t))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1]
  J = det(_jacobian_matrix)

  vx = grid_velocities(m, i, t)
  ξt = -(vx * ξx)

  ξ̂ = Metric1D(ξx * J, ξt * J)

  return (; ξ̂, J)
end

# ------------------------------------------------------------------
# Jacobian related functions
# ------------------------------------------------------------------

@inline function jacobian_matrix(m::CurvilinearGrid1D, i::Real, t::Real=0)
  return checkeps(SMatrix{1,1}(m.∂x∂ξ(i - m.nhalo, t)))
end

@inline function jacobian(m::CurvilinearGrid1D, i::Real, t::Real=0)
  return abs(m.∂x∂ξ(i - m.nhalo, t))
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::CurvilinearGrid1D, i, t::Real=0) = 0.0
# @inline centroid_velocities(m::CurvilinearGrid1D, i, t) = 0.0
# @inline node_velocities(m::CurvilinearGrid1D, i, t) = 0.0

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

function _node_coordinates!(
  coordinates::StructArray{T,1}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the node coordinates
  @inbounds for idx in domain
    cell_idx = @. idx.I - nhalo
    coordinates.x[idx] = coordinate_functions.x(cell_idx...)
  end

  return nothing
end

function _centroid_coordinates!(
  centroids::StructArray{T,1}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the centroid coordinates
  @inbounds for idx in domain
    cell_idx = @. idx.I - nhalo + 0.5
    centroids.x[idx] = coordinate_functions.x(cell_idx...)
  end

  return nothing
end
