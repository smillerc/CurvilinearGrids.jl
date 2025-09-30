
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
  discretization_scheme_name::Symbol
end

struct UniformGrid1D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid1D
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
  discretization_scheme_name::Symbol
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
  discretization_scheme_name::Symbol
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
  discretization_scheme_name::Symbol
end

"""
    CurvilinearGrid1D(x::AbstractVector{T}, discretization_scheme::Symbol; backend=CPU(), is_static=false, empty_metrics=false) where {T}

Construct a curvilinear grid in 1D using a vector of x coordinate points.
"""
function CurvilinearGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=KernelAbstractions.CPU(),
  is_static=false,
  empty_metrics=false,
  halo_coords_included=false,
) where {T}

  #
  MetricDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_metric_disc_scheme(
    discretization_scheme
  )

  limits, iterators = get_iterators(size(x), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = MetricDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes = _grid_constructor(
    x, iterators; backend=backend
  )

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  m = CurvilinearGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    iterators,
    discr_scheme,
    is_static,
    scheme_name,
  )

  if !empty_metrics
    update!(m; force=true)
  end
  return m
end

function CurvilinearGrid1D(
  (x0, x1)::NTuple, ni::Int, discretization_scheme::Symbol; kwargs...
)
  CurvilinearGrid1D(
    collect(range(x0, x1; length=ni + 1)), discretization_scheme::Symbol; kwargs...
  )
end

"""
    UniformGrid1D((x0, x1), ncells, discretization_scheme::Symbol; backend=CPU(), T=Float64, empty_metrics=false)

TBW
"""
function UniformGrid1D(
  (x0, x1),
  ncells,
  discretization_scheme::Symbol;
  backend=KernelAbstractions.CPU(),
  T=Float64,
  empty_metrics=false,
  halo_coords_included=false,
)
  ni = ncells + 1
  x = collect(T, range(x0, x1; length=ni))

  m = UniformGrid1D(
    x,
    discretization_scheme;
    backend=backend,
    empty_metrics=empty_metrics,
    halo_coords_included=halo_coords_included,
  )

  if !empty_metrics
    update!(m; force=true)
  end
  return m
end

"""
    UniformGrid1D(x::AbstractVector{T}, discretization_scheme::Symbol; backend=CPU(), is_static=true, empty_metrics=false) where {T}

TBW
"""
function UniformGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=KernelAbstractions.CPU(),
  is_static=true,
  empty_metrics=false,
  halo_coords_included=false,
) where {T}

  #
  ni = length(x)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  MetricDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_metric_disc_scheme(
    discretization_scheme
  )

  limits, iterators = get_iterators(size(x), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = MetricDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes = _grid_constructor(
    x, iterators; backend=backend
  )

  cell_center_metrics, edge_metrics = get_metric_soa_uniform1d(celldims, backend, T)

  m = UniformGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    iterators,
    discr_scheme,
    is_static,
    scheme_name,
  )

  update!(m; force=true)

  return m
end

function _grid_constructor(
  x::AbstractVector{T}, domain_iterators; backend=KernelAbstractions.CPU(), kwargs...
) where {T}
  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  centroids = StructArray((x=KernelAbstractions.zeros(backend, T, celldims),))
  coords = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
  end

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  nnodes = length(x)
  return (coords, centroids, node_velocities, nnodes)
end

"""
    SphericalGrid1D(x::AbstractVector, discretization_scheme::Symbol, snap_to_axis::Bool; backend=CPU(), is_static=false)

Construct a curvilinear grid in 1D with spherical symmetry using a vector of x coordinate points.
"""
function SphericalGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol,
  snap_to_axis::Bool;
  backend=CPU(),
  is_static=false,
  halo_coords_included=false,
  kwargs...,
) where {T}

  #
  MetricDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_metric_disc_scheme(
    discretization_scheme
  )

  limits, iterators = get_iterators(size(x), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = MetricDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes = _grid_constructor(
    x, iterators; backend=backend
  )

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  m = SphericalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    iterators,
    discr_scheme,
    snap_to_axis,
    is_static,
    scheme_name,
  )

  update!(m; force=true)
  return m
end

"""
    CylindricalGrid1D(x::AbstractVector, discretization_scheme::Symbol, snap_to_axis::Bool; backend=CPU(), is_static=false)

Construct a curvilinear grid in 1D with cylindrical symmetry using a vector of x coordinate points.
"""
function CylindricalGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol,
  snap_to_axis::Bool;
  backend=CPU(),
  is_static=false,
  halo_coords_included=false,
  kwargs...,
) where {T}

  #
  MetricDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_metric_disc_scheme(
    discretization_scheme
  )

  limits, iterators = get_iterators(size(x), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = MetricDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes = _grid_constructor(
    x, iterators; backend=backend
  )

  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  m = CylindricalGrid1D(
    coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    iterators,
    discr_scheme,
    snap_to_axis,
    is_static,
    scheme_name,
  )

  update!(m; force=true)
  return m
end

"""Update metrics after grid coordinates change"""
function update!(
  mesh::AbstractCurvilinearGrid1D; force=false, include_halo_region::Bool=false
)
  if include_halo_region
    metric_domain = mesh.iterators.cell.full
  else
    metric_domain = mesh.iterators.cell.domain
  end

  if !mesh.is_static || force
    backend = KernelAbstractions.get_backend(mesh.centroid_coordinates.x)
    _centroid_coordinates_kernel!(backend)(
      mesh.centroid_coordinates,
      mesh.node_coordinates,
      mesh.iterators.cell.domain;
      ndrange=size(mesh.iterators.cell.domain),
    )
    update_metrics!(mesh; include_halo_region=include_halo_region)
    _check_valid_metrics(mesh)
  end

  return nothing
end

function jacobian_matrix(mesh::AbstractCurvilinearGrid1D, (i,)::NTuple{1,Int})
  xξ = mesh.cell_center_metrics.x₁.ξ
  return @SMatrix [xξ[i]]
end

"""
    forward_cell_metrics(mesh::AbstractCurvilinearGrid1D, idx)

Get the forward cell metrics for a given cell index `idx`
"""
function forward_cell_metrics(mesh::AbstractCurvilinearGrid1D, idx)
  (; x=mesh.cell_center_metrics.x₁[idx],)
end

"""
    inverse_cell_metrics(mesh::AbstractCurvilinearGrid2D, idx)

Get the inverse cell metrics for a given cell index `idx`
"""
function inverse_cell_metrics(mesh::AbstractCurvilinearGrid1D, idx)
  (; ξ=mesh.cell_center_metrics.ξ[idx], ξ̂=mesh.cell_center_metrics.ξ̂[idx])
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::AbstractCurvilinearGrid1D, i, t::Real=0) = 0.0
# @inline centroid_velocities(m::AbstractCurvilinearGrid1D, i, t) = 0.0
# @inline node_velocities(m::AbstractCurvilinearGrid1D, i, t) = 0.0

# ------------------------------------------------------------------
# Coordinate Functions
# ------------------------------------------------------------------

@kernel inbounds = true function _centroid_coordinates_kernel!(
  centroids::StructArray{T,1}, coords::StructArray{T,1}, domain
) where {T}
  idx = @index(Global)
  didx = domain[idx]
  i, = didx.I

  x = coords.x

  centroids.x[didx] = 0.5(x[i] + x[i + 1])
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

    edge_metrics_valid =
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain]))
  end

  if !edge_metrics_valid
    error("Invalid edge metrics found")
  end

  if !centroid_metrics_valid
    error("Invalid centroid metrics found")
  end

  return nothing
end
