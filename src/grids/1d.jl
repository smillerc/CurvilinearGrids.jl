
"""
CurvilinearGrid1D

# Fields
 - `x`: Node function, e.g. x -> f(ξ)
 - `∂x∂ξ`: Derivative of x wrt ξ; ∂x∂ξ(ξ)
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CurvilinearGrid1D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  is_static::Bool
  is_orthogonal::Bool
end

struct SphericalGrid1D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  snap_to_axis::Bool
  is_static::Bool
  is_orthogonal::Bool
end

struct CylindricalGrid1D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  snap_to_axis::Bool
  is_static::Bool
  is_orthogonal::Bool
end

"""
    CurvilinearGrid1D(x::Function, (n_ξ,), nhalo)

Create a `CurvilinearGrid1D` with a function `x(ξ)` and `nhalo` halo cells. `n_ξ` is the
total number of nodes/vertices (not including halo).
"""
function CurvilinearGrid1D(
  x::AbstractVector{T}, nhalo; backend=CPU(), discretization_scheme=:MEG6, is_static=false
) where {T}
  ni = length(x)
  ncells = ni - 1
  ilo = nhalo + 1
  limits = (node=(ilo=ilo, ihi=ni + nhalo), cell=(ilo=ilo, ihi=ncells + nhalo))

  nodeCI = CartesianIndices((ni + 2nhalo,))
  cellCI = CartesianIndices((ncells + 2nhalo,))

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  centroids = StructArray((x=KernelAbstractions.zeros(backend, T, celldims),))
  coords = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
  end

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

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

  m = CurvilinearGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    ni,
    limits,
    domain_iterators,
    discr_scheme,
    is_static,
    true,
  )

  update!(m; force=true)
  return m
end

function SphericalGrid1D(
  x::AbstractVector{T},
  nhalo::Int,
  snap_to_axis::Bool;
  backend=CPU(),
  discretization_scheme=:MEG6,
  is_static=false,
) where {T}
  ni = length(x)
  ncells = ni - 1
  ilo = nhalo + 1
  limits = (node=(ilo=ilo, ihi=ni + nhalo), cell=(ilo=ilo, ihi=ncells + nhalo))

  nodeCI = CartesianIndices((ni + 2nhalo,))
  cellCI = CartesianIndices((ncells + 2nhalo,))

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  centroids = StructArray((x=KernelAbstractions.zeros(backend, T, celldims),))
  coords = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
  end

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

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

  m = SphericalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    ni,
    limits,
    domain_iterators,
    discr_scheme,
    snap_to_axis,
    is_static,
    true,
  )

  update!(m; force=true)
  return m
end

function CylindricalGrid1D(
  x::AbstractVector{T},
  nhalo::Int,
  snap_to_axis::Bool,
  ;
  backend=CPU(),
  discretization_scheme=:MEG6,
  is_static=false,
) where {T}
  ni = length(x)
  ncells = ni - 1
  ilo = nhalo + 1
  limits = (node=(ilo=ilo, ihi=ni + nhalo), cell=(ilo=ilo, ihi=ncells + nhalo))

  nodeCI = CartesianIndices((ni + 2nhalo,))
  cellCI = CartesianIndices((ncells + 2nhalo,))

  domain_iterators = get_node_cell_iterators(nodeCI, cellCI, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  centroids = StructArray((x=KernelAbstractions.zeros(backend, T, celldims),))
  coords = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
  end

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

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

  m = CylindricalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    ni,
    limits,
    domain_iterators,
    discr_scheme,
    snap_to_axis,
    is_static,
    true,
  )

  update!(m; force=true)
  return m
end

"""Update metrics after grid coordinates change"""
function update!(mesh::CurvilinearGrid1D; force=false)
  if !mesh.is_static || force
    _centroid_coordinates!(
      mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
    )
    update_metrics!(mesh)
    _check_valid_metrics(mesh)
  end

  return nothing
end

function update_metrics!(mesh::AbstractCurvilinearGrid1D, t::Real=0)
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

# # ------------------------------------------------------------------
# # Grid Metrics
# # ------------------------------------------------------------------

# @inline function metrics(m::AbstractCurvilinearGrid1D, i::Real, t::Real=0)
#   _jacobian_matrix = checkeps(m._∂x∂ξ(i - m.nhalo, t))
#   inv_jacobian_matrix = inv(_jacobian_matrix)
#   ξx = inv_jacobian_matrix[1]
#   J = det(_jacobian_matrix)

#   vx = grid_velocities(m, i, t)
#   ξt = -(vx * ξx)

#   ξ = Metric1D(ξx, ξt)

#   return (; ξ, J)
# end

# # ------------------------------------------------------------------
# # Conservative Grid Metrics; e.g. ξ̂x = ξx * J
# # ------------------------------------------------------------------

# @inline function conservative_metrics(m::AbstractCurvilinearGrid1D, i::Real, t::Real=0)
#   _jacobian_matrix = checkeps(m._∂x∂ξ(i - m.nhalo, t))
#   inv_jacobian_matrix = inv(_jacobian_matrix)
#   ξx = inv_jacobian_matrix[1]
#   J = det(_jacobian_matrix)

#   vx = grid_velocities(m, i, t)
#   ξt = -(vx * ξx)

#   ξ̂ = Metric1D(ξx * J, ξt * J)

#   return (; ξ̂, J)
# end

# # ------------------------------------------------------------------
# # Jacobian related functions
# # ------------------------------------------------------------------

# @inline function jacobian_matrix(m::AbstractCurvilinearGrid1D, i::Real, t::Real=0)
#   return checkeps(SMatrix{1,1}(m.∂x∂ξ(i - m.nhalo, t)))
# end

# @inline function jacobian(m::AbstractCurvilinearGrid1D, i::Real, t::Real=0)
#   return abs(m.∂x∂ξ(i - m.nhalo, t))
# end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::AbstractCurvilinearGrid1D, i, t::Real=0) = 0.0
# @inline centroid_velocities(m::AbstractCurvilinearGrid1D, i, t) = 0.0
# @inline node_velocities(m::AbstractCurvilinearGrid1D, i, t) = 0.0

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

function _centroid_coordinates!(
  centroids::StructArray{T,1}, coords::StructArray{T,1}, domain
) where {T}
  x = coords.x
  # Populate the centroid coordinates
  for idx in domain
    i, = idx.I
    centroids.x[idx] = 0.5(x[i] + x[i + 1])
  end

  return nothing
end

function _check_valid_metrics(mesh::AbstractCurvilinearGrid1D)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, 1, -1)

  @views begin
    centroid_metrics_valid =
      all(isfinite.(mesh.cell_center_metrics.J[domain])) &&
      all(mesh.cell_center_metrics.J[domain] .> 0) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.ξ[domain]))

    # edge_metrics_valid =
    # all(isfinite.(mesh.edge_metrics.i₊½.J[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
    #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain])) &&

  end

  # if !edge_metrics_valid
  #   error("Invalid edge metrics found")
  # end

  if !centroid_metrics_valid
    error("Invalid centroid metrics found")
  end

  return nothing
end
