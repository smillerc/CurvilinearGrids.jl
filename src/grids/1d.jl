
struct CurvilinearGrid1D{CO,CE,NV,EM,CM,DL,CI,TI,DS} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  tiles::TI
  discretization_scheme::DS
  onbc::@NamedTuple{ilo::Bool, ihi::Bool}
  is_static::Bool
  is_orthogonal::Bool
end

struct SphericalGrid1D{CO,CE,NV,EM,CM,DL,CI,TI,DS} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  tiles::TI
  discretization_scheme::DS
  onbc::@NamedTuple{ilo::Bool, ihi::Bool}
  snap_to_axis::Bool
  is_static::Bool
  is_orthogonal::Bool
end

struct CylindricalGrid1D{CO,CE,NV,EM,CM,DL,CI,TI,DS} <: AbstractCurvilinearGrid1D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::Int
  domain_limits::DL
  iterators::CI
  tiles::TI
  discretization_scheme::DS
  onbc::@NamedTuple{ilo::Bool, ihi::Bool}
  snap_to_axis::Bool
  is_static::Bool
  is_orthogonal::Bool
end

"""
    CurvilinearGrid1D(x, nhalo; backend=CPU(), discretization_scheme=:MEG6, on_bc=nothing, is_static=false, tiles=nothing)

Construct a curvilinear grid in 1D using a vector of x coordinate points.
"""
function CurvilinearGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  tiles=nothing,
) where {T}

  #
  if Symbol(uppercase("$discretization_scheme")) === :MEG6 ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
    nhalo = 5
  else
    error("Only MontoneExplicitGradientScheme6thOrder or MEG6 is supported for now")
  end

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

  discr_scheme = MetricDiscretizationScheme(;
    use_cache=true, celldims=size(domain_iterators.cell.full), backend=backend, T=T
  )

  if isnothing(on_bc)
    _on_bc = (ilo=true, ihi=true)
  else
    _on_bc = on_bc
  end

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
    tiles,
    discr_scheme,
    _on_bc,
    is_static,
    true,
  )

  update!(m; force=true)
  return m
end

"""
    SphericalGrid1D(x, nhalo, snap_to_axis::Bool; backend=CPU(), discretization_scheme=:MEG6, on_bc=nothing, is_static=false, tiles=nothing)

Construct a curvilinear grid in 1D with spherical symmetry using a vector of x coordinate points.
"""
function SphericalGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol,
  snap_to_axis::Bool;
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  tiles=nothing,
) where {T}

  #
  if Symbol(uppercase("$discretization_scheme")) === :MEG6 ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
    nhalo = 5
  else
    error("Only MontoneExplicitGradientScheme6thOrder or MEG6 is supported for now")
  end

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

  discr_scheme = MetricDiscretizationScheme(;
    use_cache=true, celldims=size(domain_iterators.cell.full), backend=backend, T=T
  )

  if isnothing(on_bc)
    _on_bc = (ilo=true, ihi=true)
  else
    _on_bc = on_bc
  end

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
    tiles,
    discr_scheme,
    _on_bc,
    snap_to_axis,
    is_static,
    true,
  )

  update!(m; force=true)
  return m
end

"""
    CylindricalGrid1D(x, nhalo, snap_to_axis::Bool; backend=CPU(), discretization_scheme=:MEG6, on_bc=nothing, is_static=false, tiles=nothing)

Construct a curvilinear grid in 1D with cylindrical symmetry using a vector of x coordinate points.
"""
function CylindricalGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol,
  snap_to_axis::Bool;
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  tiles=nothing,
) where {T}

  #
  if Symbol(uppercase("$discretization_scheme")) === :MEG6 ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
    nhalo = 5
  else
    error("Only MontoneExplicitGradientScheme6thOrder or MEG6 is supported for now")
  end

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

  discr_scheme = MetricDiscretizationScheme(;
    use_cache=true, celldims=size(domain_iterators.cell.full), backend=backend, T=T
  )

  if isnothing(on_bc)
    _on_bc = (ilo=true, ihi=true)
  else
    _on_bc = on_bc
  end

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
    tiles,
    discr_scheme,
    _on_bc,
    snap_to_axis,
    is_static,
    true,
  )

  update!(m; force=true)
  return m
end

"""Update metrics after grid coordinates change"""
function update!(mesh::AbstractCurvilinearGrid1D; force=false)
  if !mesh.is_static || force
    _centroid_coordinates!(
      mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
    )
    update_metrics!(mesh)
    _check_valid_metrics(mesh)
  end

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
  # domain = mesh.iterators.cell.domain
  # i₊½_domain = expand(domain, 1, -1)

  # @views begin
  #   centroid_metrics_valid =
  #     all(isfinite.(mesh.cell_center_metrics.J[domain])) &&
  #     all(mesh.cell_center_metrics.J[domain] .> 0) &&
  #     all(isfinite.(mesh.cell_center_metrics.ξ.x₁[domain])) &&
  #     all(isfinite.(mesh.cell_center_metrics.x₁.ξ[domain]))

  #   # edge_metrics_valid =
  #   # all(isfinite.(mesh.edge_metrics.i₊½.J[i₊½_domain])) &&
  #   #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
  #   #     all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain])) &&

  # end

  # # if !edge_metrics_valid
  # #   error("Invalid edge metrics found")
  # # end

  # if !centroid_metrics_valid
  #   error("Invalid centroid metrics found")
  # end

  return nothing
end
