
"""
CylindricalGrid2D

# Fields
 - `x`: Node function; e.g., x(i,j)
 - `y`: Node function; e.g., y(i,j)
 - `jacobian_matrix_func`: jacobian matrix, e.g., J(i,j)
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CylindricalGrid2D{CO,CE,EMC,NV,EM,CM,CI,CF,JF} <: AbstractCurvilinearGrid
  node_coordinates::CO
  centroid_coordinates::CE
  edge_midpoint_coordinates::EMC
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{2,Int}
  iterators::CI
  snap_to_axis::Bool
  _coordinate_funcs::CF
  _jacobian_matrix_func::JF
end

"""
    CylindricalGrid2D(r, z, (ni, nj), nhalo, snap_to_axis; T=Float64, backend=CPU())

"""
function CylindricalGrid2D(
  r::Function, z::Function, (ni, nj), nhalo, snap_to_axis; T=Float64, backend=CPU()
)
  dim = 2
  check_nargs(r, dim, :r)
  check_nargs(z, dim, :z)
  test_coord_func(r, dim, :r)
  test_coord_func(z, dim, :z)

  rz(i, j) = @SVector [r(i, j), z(i, j)]
  function jacobian_matrix_func(i, j, t)
    return ForwardDiff.jacobian(x -> rz(x[1], x[2]), @SVector [i, j])
  end

  nnodes = (ni, nj)
  ncells = nnodes .- 1

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)
  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  coordinate_funcs = (; r, z)
  centroids = StructArray((
    r=KernelAbstractions.zeros(backend, T, celldims),
    z=KernelAbstractions.zeros(backend, T, celldims),
  ))

  coords = StructArray((
    r=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  edge_coords = (
    i₊½=StructArray((
      r=KernelAbstractions.zeros(backend, T, celldims),
      z=KernelAbstractions.zeros(backend, T, celldims),
    )),
    j₊½=StructArray((
      r=KernelAbstractions.zeros(backend, T, celldims),
      z=KernelAbstractions.zeros(backend, T, celldims),
    )),
  )

  node_velocities = StructArray((
    r=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  m = CylindricalGrid2D(
    coords,
    centroids,
    edge_coords,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    domain_iterators,
    snap_to_axis,
    coordinate_funcs,
    jacobian_matrix_func,
  )

  update_coordinates(m)
  update_metrics!(m)
  check_for_invalid_metrics(m)
  return m
end

function update_coordinates(mesh::CylindricalGrid2D)
  update_coordinates(
    mesh.centroid_coordinates,
    mesh.node_coordinates,
    mesh.edge_midpoint_coordinates,
    mesh._coordinate_funcs,
    mesh.iterators,
    mesh.nhalo,
    mesh.snap_to_axis,
  )

  return nothing
end

function update_coordinates(
  centroids, coords, edge_coords, coordinate_funcs, domain_iterators, nhalo, snap_to_axis
)
  _rz_centroid_coordinates!(centroids, coordinate_funcs, domain_iterators.cell.full, nhalo)
  _rz_node_coordinates!(coords, coordinate_funcs, domain_iterators.node.full, nhalo)

  if snap_to_axis
    _check_nodes_along_axis(coords, domain_iterators.node.domain)
  end

  _rz_edge_coordinates!(edge_coords, coordinate_funcs, domain_iterators.cell.full, nhalo)

  return nothing
end

function update_metrics!(mesh::CylindricalGrid2D, t::Real=0)

  # cell metrics
  @inbounds for idx in mesh.iterators.cell.full
    cell_idx = idx.I .+ 0.5
    # @unpack J, ξ, η, x, y = metrics(mesh, cell_idx, t)
    @unpack J, ξ, η = metrics(mesh, cell_idx, t)

    mesh.cell_center_metrics.ξ.x₁[idx] = ξ.x₁
    mesh.cell_center_metrics.ξ.x₂[idx] = ξ.x₂
    mesh.cell_center_metrics.ξ.t[idx] = ξ.t
    mesh.cell_center_metrics.η.x₁[idx] = η.x₁
    mesh.cell_center_metrics.η.x₂[idx] = η.x₂
    mesh.cell_center_metrics.η.t[idx] = η.t

    # mesh.cell_center_inv_metrics.xξ[idx] = x.ξ
    # mesh.cell_center_inv_metrics.yξ[idx] = y.ξ
    # mesh.cell_center_inv_metrics.xη[idx] = x.η
    # mesh.cell_center_inv_metrics.yη[idx] = y.η

    mesh.cell_center_metrics.J[idx] = J
  end

  # i₊½ conserved metrics
  @inbounds for idx in mesh.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i₊½, j)
    @unpack ξ̂, η̂, J = conservative_metrics(mesh, (i + 1 / 2, j), t)

    mesh.edge_metrics.i₊½.ξ̂.x₁[idx] = ξ̂.x₁
    mesh.edge_metrics.i₊½.ξ̂.x₂[idx] = ξ̂.x₂
    mesh.edge_metrics.i₊½.ξ̂.t[idx] = ξ̂.t
    mesh.edge_metrics.i₊½.η̂.x₁[idx] = η̂.x₁
    mesh.edge_metrics.i₊½.η̂.x₂[idx] = η̂.x₂
    mesh.edge_metrics.i₊½.η̂.t[idx] = η̂.t
    mesh.edge_metrics.i₊½.J[idx] = J
  end

  # j₊½ conserved metrics
  @inbounds for idx in mesh.iterators.cell.full
    i, j = idx.I .+ 0.5 # centroid index

    # get the conserved metrics at (i, j₊½)
    @unpack ξ̂, η̂, J = conservative_metrics(mesh, (i, j + 1 / 2), t)

    mesh.edge_metrics.j₊½.ξ̂.x₁[idx] = ξ̂.x₁
    mesh.edge_metrics.j₊½.ξ̂.x₂[idx] = ξ̂.x₂
    mesh.edge_metrics.j₊½.ξ̂.t[idx] = ξ̂.t
    mesh.edge_metrics.j₊½.η̂.x₁[idx] = η̂.x₁
    mesh.edge_metrics.j₊½.η̂.x₂[idx] = η̂.x₂
    mesh.edge_metrics.j₊½.η̂.t[idx] = η̂.t
    mesh.edge_metrics.j₊½.J[idx] = J
  end

  return nothing
end

# # ------------------------------------------------------------------
# # Grid Metrics
# # ------------------------------------------------------------------

@inline function metrics(mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  _jacobian_matrix = checkeps(mesh._jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo, t))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξr = inv_jacobian_matrix[1, 1]
  ξz = inv_jacobian_matrix[1, 2]
  ηr = inv_jacobian_matrix[2, 1]
  ηz = inv_jacobian_matrix[2, 2]

  J = det(_jacobian_matrix)

  vr, vz = grid_velocities(mesh, (i, j), t)
  ξt = -(vr * ξr + vz * ξz)
  ηt = -(vr * ηr + vz * ηz)

  ξ = Metric2D(ξr, ξz, ξt)
  η = Metric2D(ηr, ηz, ηt)

  return (; ξ, η, J)
end

# # ------------------------------------------------------------------
# # Conservative Grid Metrics; e.g. ξ̂x = ξr * J
# # ------------------------------------------------------------------

@inline function conservative_metrics(
  mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t::Real=0
)
  _jacobian_matrix = checkeps(mesh._jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo, t))
  inv_jacobian_matrix = inv(_jacobian_matrix)
  ξr = inv_jacobian_matrix[1, 1]
  ξz = inv_jacobian_matrix[1, 2]
  ηr = inv_jacobian_matrix[2, 1]
  ηz = inv_jacobian_matrix[2, 2]

  J = det(_jacobian_matrix)

  vr, vz = grid_velocities(mesh, (i, j), t)
  ξt = -(vr * ξr + vz * ξz)
  ηt = -(vr * ηr + vz * ηz)

  ξ̂ = Metric2D(ξr * J, ξz * J, ξt * J)
  η̂ = Metric2D(ηr * J, ηz * J, ηt * J)

  return (; ξ̂, η̂, J)
end

# # ------------------------------------------------------------------
# # Jacobian related functions
# # ------------------------------------------------------------------
function jacobian_matrix(mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  return checkeps(mesh._jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo, t))
end

function jacobian(mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  return det(jacobian_matrix(mesh, (i, j), t))
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t::Real=0) = (0.0, 0.0)
# @inline centroid_velocities(mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline node_velocities(mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

function _rz_node_coordinates!(
  coordinates::StructArray{T,2}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the node coordinates
  for idx in domain
    cell_idx = @. idx.I - nhalo
    coordinates.r[idx] = coordinate_functions.r(cell_idx...)
    coordinates.z[idx] = coordinate_functions.z(cell_idx...)
  end

  return nothing
end

function _rz_centroid_coordinates!(
  centroids::StructArray{T,2}, coordinate_functions, domain, nhalo
) where {T}

  # Populate the centroid coordinates
  for idx in domain
    cell_idx = @. idx.I - nhalo + 0.5
    centroids.r[idx] = coordinate_functions.r(cell_idx...)
    centroids.z[idx] = coordinate_functions.z(cell_idx...)
  end

  return nothing
end

function _rz_edge_coordinates!(edge_coords, coordinate_functions, domain, nhalo)

  # Populate the centroid coordinates
  for idx in domain
    cell_idx = @. idx.I - nhalo + 0.5

    i₊½ = cell_idx .+ (0.5, 0)
    j₊½ = cell_idx .+ (0, 0.5)
    edge_coords.i₊½.r[idx] = coordinate_functions.r(i₊½...)
    edge_coords.i₊½.z[idx] = coordinate_functions.z(i₊½...)

    edge_coords.j₊½.r[idx] = coordinate_functions.r(j₊½...)
    edge_coords.j₊½.z[idx] = coordinate_functions.z(j₊½...)
  end

  return nothing
end

function _check_nodes_along_axis(nodes, inner_domain)
  @views begin
    axis_domain = inner_domain[1, :]
    if any(!iszero(nodes.r[axis_domain]))
      error(
        "Nodes not aligned to axis of symmetry (`snap_to_axis = true`). Set snap_to_axis = false to disable checks",
      )
    end
  end
end
