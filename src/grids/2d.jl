
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

"""
    CurvilinearGrid2D(x, y, nhalo::Int, discretization_scheme=:MEG6; backend=CPU())

"""
function CurvilinearGrid2D(
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  nhalo::Int,
  discretization_scheme=:MEG6;
  backend=CPU(),
) where {T}

  #
  @assert size(x) == size(y)

  nnodes = size(x)
  ni, nj = nnodes

  ncells = nnodes .- 1
  ni_cells, nj_cells = ncells
  lo = nhalo + 1
  limits = (
    node=(ilo=lo, ihi=ni + nhalo, jlo=lo, jhi=nj + nhalo),
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

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
    copy!(coords.y[domain_iterators.node.domain], y)
  end

  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
  ))

  _centroid_coordinates!(centroids, coords, domain_iterators.cell.domain)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  # if discretization_scheme === :MEG6 || discretization_scheme === :MonotoneExplicit6thOrder
  #   if nhalo < 4
  #     error("`nhalo` must = 4 when using the MEG6 discretization scheme")
  #   end

  discr_scheme = MetricDiscretizationSchemes.MonotoneExplicit6thOrderDiscretization(
    domain_iterators.cell.full
  )

  # else
  #   error("Unknown discretization scheme to compute the conserved metrics")
  # end

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
  check_valid_metrics(mesh)
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

function update!(mesh::CurvilinearGrid2D)
  _centroid_coordinates!(
    mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
  )
  update_metrics!(mesh)
  check_valid_metrics(mesh)
  return nothing
end

function update_metrics!(mesh::CurvilinearGrid2D, t::Real=0)
  # Update the metrics within the non-halo region, e.g., the domain
  domain = mesh.iterators.cell.domain

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

function check_valid_metrics(mesh::CurvilinearGrid2D)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, 1, -1)
  j₊½_domain = expand(domain, 2, -1)

  @views begin
    centroid_metrics_valid = all(isfinite.(mesh.cell_center_metrics.J[domain])) # &&
    #     all(isfinite.(mesh.cell_center_metrics.ξ.x₁[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.ξ.x₂[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.η.x₁[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.η.x₂[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₁.ξ[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₁.η[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₂.ξ[domain])) &&
    #     all(isfinite.(mesh.cell_center_metrics.x₂.η[domain]))

    # edge_metrics_valid =
    # all(isfinite.(mesh.edge_metrics.i₊½.J[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₂[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₁[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₂[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.η̂.t[i₊½_domain])) &&
    # all(isfinite.(mesh.edge_metrics.j₊½.J[j₊½_domain])) # &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₁[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₂[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.t[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₁[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₂[j₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.j₊½.η̂.t[j₊½_domain]))
  end

  # if !edge_metrics_valid
  #   error("Invalid edge metrics found")
  # end

  if !centroid_metrics_valid
    error("Invalid centroid metrics found")
  end

  return nothing
end