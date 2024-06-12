
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
struct CurvilinearGrid2D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid2D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{2,Int}
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  # _coordinate_funcs::CF
  # _jacobian_matrix_func::JF
end

function CurvilinearGrid2D(
  x::AbstractArray{T,2}, y::AbstractArray{T,2}, nhalo::Int; backend=CPU()
) where {T}
  nnodes = size(x) .- 2nhalo
  n_ξ, n_η = nnodes

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

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))
  copy!(coords.x, x)
  copy!(coords.y, y)

  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
  ))

  _centroid_coordinates!(centroids, coords, domain_iterators.cell.domain)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  discr_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
    domain_iterators.cell.full
  )

  mesh = CurvilinearGrid2D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    domain_iterators,
    discr_scheme,
    # coordinate_funcs,
    # jacobian_matrix_func,
  )

  update_metrics!(mesh)
  # check_for_invalid_metrics(mesh)
  return mesh
end

# function CurvilinearGrid2D(
#   x::Function, y::Function, (n_ξ, n_η), nhalo; T=Float64, backend=CPU()
# )
#   dim = 2
#   check_nargs(x, dim, :x)
#   check_nargs(y, dim, :y)
#   test_coord_func(x, dim, :x)
#   test_coord_func(y, dim, :y)

#   xy(i, j) = @SVector [x(i, j), y(i, j)]
#   function jacobian_matrix_func(i, j, t)
#     return ForwardDiff.jacobian(x -> xy(x[1], x[2]), @SVector [i, j])
#   end
#   # jacobian_matrix_func = _setup_jacobian_func(x, y)
#   nnodes = (n_ξ, n_η)
#   ncells = nnodes .- 1
#   ni_cells, nj_cells = ncells
#   lo = nhalo + 1
#   limits = (
#     node=(ilo=lo, ihi=n_ξ + nhalo, jlo=lo, jhi=n_η + nhalo),
#     cell=(ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo),
#   )

#   nodeCI = CartesianIndices(nnodes .+ 2nhalo)
#   cellCI = CartesianIndices(ncells .+ 2nhalo)

#   domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)
#   celldims = size(domain_iterators.cell.full)
#   nodedims = size(domain_iterators.node.full)

#   cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

#   coordinate_funcs = (; x, y)
#   centroids = StructArray((
#     x=KernelAbstractions.zeros(backend, T, celldims),
#     y=KernelAbstractions.zeros(backend, T, celldims),
#   ))
#   _centroid_coordinates!(centroids, coordinate_funcs, domain_iterators.cell.full, nhalo)

#   coords = StructArray((
#     x=KernelAbstractions.zeros(backend, T, nodedims),
#     y=KernelAbstractions.zeros(backend, T, nodedims),
#   ))
#   _node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

#   node_velocities = StructArray((
#     x=KernelAbstractions.zeros(backend, T, nodedims),
#     y=KernelAbstractions.zeros(backend, T, nodedims),
#   ))

#   discr_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
#     domain_iterators.cell.full
#   )

#   mesh = CurvilinearGrid2D(
#     coords,
#     centroids,
#     node_velocities,
#     edge_metrics,
#     cell_center_metrics,
#     nhalo,
#     nnodes,
#     limits,
#     domain_iterators,
#     discr_scheme,
#     coordinate_funcs,
#     jacobian_matrix_func,
#   )

#   update_metrics!(mesh)
#   # check_for_invalid_metrics(mesh)
#   return mesh
# end

function update_metrics!(mesh::CurvilinearGrid2D, t::Real=0)
  # Update the metrics within the non-halo region, e.g., the domain
  domain = mesh.iterators.cell.domain

  # # cell metrics
  # @inbounds for idx in mesh.iterators.cell.full
  #   cell_idx = idx.I .+ 0.5
  #   # @unpack J, ξ, η, x, y = metrics(mesh, cell_idx, t)
  #   @unpack J, ξ, η = metrics(mesh, cell_idx, t)

  #   mesh.cell_center_metrics.ξ.x₁[idx] = ξ.x₁
  #   mesh.cell_center_metrics.ξ.x₂[idx] = ξ.x₂
  #   mesh.cell_center_metrics.ξ.t[idx] = ξ.t
  #   mesh.cell_center_metrics.η.x₁[idx] = η.x₁
  #   mesh.cell_center_metrics.η.x₂[idx] = η.x₂
  #   mesh.cell_center_metrics.η.t[idx] = η.t

  #   # mesh.cell_center_inv_metrics.xξ[idx] = x.ξ
  #   # mesh.cell_center_inv_metrics.yξ[idx] = y.ξ
  #   # mesh.cell_center_inv_metrics.xη[idx] = x.η
  #   # mesh.cell_center_inv_metrics.yη[idx] = y.η

  #   mesh.cell_center_metrics.J[idx] = J
  # end

  # # i₊½ conserved metrics
  # @inbounds for idx in mesh.iterators.cell.full
  #   i, j = idx.I .+ 0.5 # centroid index

  #   # get the conserved metrics at (i₊½, j)
  #   @unpack ξ̂, η̂, J = conservative_metrics(mesh, (i + 1 / 2, j), t)

  #   mesh.edge_metrics.i₊½.ξ̂.x₁[idx] = ξ̂.x₁
  #   mesh.edge_metrics.i₊½.ξ̂.x₂[idx] = ξ̂.x₂
  #   mesh.edge_metrics.i₊½.ξ̂.t[idx] = ξ̂.t
  #   mesh.edge_metrics.i₊½.η̂.x₁[idx] = η̂.x₁
  #   mesh.edge_metrics.i₊½.η̂.x₂[idx] = η̂.x₂
  #   mesh.edge_metrics.i₊½.η̂.t[idx] = η̂.t
  #   mesh.edge_metrics.i₊½.J[idx] = J
  # end

  # # j₊½ conserved metrics
  # @inbounds for idx in mesh.iterators.cell.full
  #   i, j = idx.I .+ 0.5 # centroid index

  #   # get the conserved metrics at (i, j₊½)
  #   @unpack ξ̂, η̂, J = conservative_metrics(mesh, (i, j + 1 / 2), t)

  #   mesh.edge_metrics.j₊½.ξ̂.x₁[idx] = ξ̂.x₁
  #   mesh.edge_metrics.j₊½.ξ̂.x₂[idx] = ξ̂.x₂
  #   mesh.edge_metrics.j₊½.ξ̂.t[idx] = ξ̂.t
  #   mesh.edge_metrics.j₊½.η̂.x₁[idx] = η̂.x₁
  #   mesh.edge_metrics.j₊½.η̂.x₂[idx] = η̂.x₂
  #   mesh.edge_metrics.j₊½.η̂.t[idx] = η̂.t
  #   mesh.edge_metrics.j₊½.J[idx] = J
  # end

  MetricDiscretizationSchemes.update_metrics!(
    mesh.discretization_scheme,
    mesh.centroid_coordinates,
    mesh.cell_center_metrics,
    mesh.edge_metrics,
    domain,
  )

  return nothing
end

# ------------------------------------------------------------------
# Grid Metrics
# ------------------------------------------------------------------

@inline function metrics(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  _jacobian_matrix = checkeps(mesh._jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo, t))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  J = det(_jacobian_matrix)

  vx, vy = grid_velocities(mesh, (i, j), t)
  ξt = -(vx * ξx + vy * ξy)
  ηt = -(vx * ηx + vy * ηy)

  ξ = Metric2D(ξx, ξy, ξt)
  η = Metric2D(ηx, ηy, ηt)

  return (; ξ, η, J)
end

# ------------------------------------------------------------------
# Conservative Grid Metrics; e.g. ξ̂x = ξx * J
# ------------------------------------------------------------------

@inline function conservative_metrics(
  mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0
)
  _jacobian_matrix = checkeps(mesh._jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo, t))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξx = inv_jacobian_matrix[1, 1]
  ξy = inv_jacobian_matrix[1, 2]
  ηx = inv_jacobian_matrix[2, 1]
  ηy = inv_jacobian_matrix[2, 2]

  J = det(_jacobian_matrix)

  vx, vy = grid_velocities(mesh, (i, j), t)
  ξt = -(vx * ξx + vy * ξy)
  ηt = -(vx * ηx + vy * ηy)

  ξ̂ = Metric2D(ξx * J, ξy * J, ξt * J)
  η̂ = Metric2D(ηx * J, ηy * J, ηt * J)

  return (; ξ̂, η̂, J)
end

# ------------------------------------------------------------------
# Jacobian related functions
# ------------------------------------------------------------------
function jacobian_matrix(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  return checkeps(mesh._jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo, t))
end

function jacobian(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  return det(jacobian_matrix(mesh, (i, j), t))
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0) = (0.0, 0.0)
# @inline centroid_velocities(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline node_velocities(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

# function _node_coordinates!(
#   coordinates::StructArray{T,2}, coordinate_functions, domain, nhalo
# ) where {T}

#   # Populate the node coordinates
#   @inbounds for idx in domain
#     cell_idx = @. idx.I - nhalo
#     coordinates.x[idx] = coordinate_functions.x(cell_idx...)
#     coordinates.y[idx] = coordinate_functions.y(cell_idx...)
#   end

#   return nothing
# end

function _centroid_coordinates!(
  centroids::StructArray{T,2}, coords::StructArray{T,2}, domain
) where {T}
  x = coords.x
  y = coords.y
  # Populate the centroid coordinates
  @inbounds for idx in domain
    i, j = idx.I
    centroids.x[idx] = 0.25(x[i, j] + x[i + 1, j] + x[i + 1, j + 1] + x[i, j + 1])
    centroids.y[idx] = 0.25(y[i, j] + y[i + 1, j] + y[i + 1, j + 1] + y[i, j + 1])
  end

  return nothing
end
