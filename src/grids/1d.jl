
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
  onbc::@NamedTuple{ilo::Bool, ihi::Bool}
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
  onbc::@NamedTuple{ilo::Bool, ihi::Bool}
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
  backend=CPU(),
  is_static=false,
  empty_metrics=false,
) where {T}
  m = CurvilinearGrid1D(
    _grid_constructor(
      x,
      :curvilinear,
      discretization_scheme;
      backend=backend,
      is_static=is_static,
      empty_metrics=empty_metrics,
    )...,
  )

  if !empty_metrics
    update!(m; force=true)
  end
  return m
end

function CurvilinearGrid1D((x0, x1), ni, discretization_scheme::Symbol; kwargs...)
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
  backend=CPU(),
  T=Float64,
  empty_metrics=false,
)
  ni = ncells + 1
  x = collect(T, range(x0, x1; length=ni))

  m = UniformGrid1D(
    _grid_constructor(
      x, :uniform, discretization_scheme; backend=backend, empty_metrics=empty_metrics
    )...,
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
  backend=CPU(),
  is_static=true,
  empty_metrics=false,
) where {T}

  #
  ni = length(x)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  m = UniformGrid1D(
    _grid_constructor(
      x,
      :uniform,
      discretization_scheme;
      backend=backend,
      is_static=is_static,
      empty_metrics=empty_metrics,
    )...,
  )

  update!(m; force=true)

  return m
end

function _grid_constructor(
  x::AbstractVector{T},
  tag::Symbol,
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=false,
  empty_metrics=false,
  kwargs...,
) where {T}

  #
  use_symmetric_conservative_metric_scheme = false
  nhalo = nhalo_lookup[scheme_name]

  scheme_name = Symbol(uppercase("$discretization_scheme"))
  if scheme_name === :MEG6 ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder

  elseif scheme_name === :MEG6_SYMMETRIC ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder

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

  if tag == :curvilinear
    if empty_metrics
      cell_center_metrics, edge_metrics = (nothing, nothing)
    else
      cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)
    end
  else
    if empty_metrics
      cell_center_metrics, edge_metrics = (nothing, nothing)
    else
      cell_center_metrics, edge_metrics = get_metric_soa_uniform1d(celldims, backend, T)
    end
  end

  centroids = StructArray((x=KernelAbstractions.zeros(backend, T, celldims),))
  coords = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  @views begin
    copy!(coords.x[domain_iterators.node.domain], x)
  end

  node_velocities = StructArray((x=KernelAbstractions.zeros(backend, T, nodedims),))

  discr_scheme = MetricDiscretizationScheme(;
    use_cache=true,
    celldims=size(domain_iterators.cell.full),
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=false,
  )

  return (
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
    scheme_name,
  )
end

"""
    SphericalGrid1D(x::AbstractVector, discretization_scheme::Symbol, snap_to_axis::Bool; backend=CPU(), on_bc=nothing, is_static=false)

Construct a curvilinear grid in 1D with spherical symmetry using a vector of x coordinate points.
"""
function SphericalGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol,
  snap_to_axis::Bool;
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  kwargs...,
) where {T}

  #

  scheme_name = Symbol(uppercase("$discretization_scheme"))
  nhalo = nhalo_lookup[scheme_name]

  if scheme_name === :MEG6 ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
  elseif scheme_name === :MEG6_SYMMETRIC ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
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
    use_cache=true,
    celldims=size(domain_iterators.cell.full),
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=false,
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
    discr_scheme,
    _on_bc,
    snap_to_axis,
    is_static,
    scheme_name,
  )

  update!(m; force=true)
  return m
end

"""
    CylindricalGrid1D(x::AbstractVector, discretization_scheme::Symbol, snap_to_axis::Bool; backend=CPU(), on_bc=nothing, is_static=false)

Construct a curvilinear grid in 1D with cylindrical symmetry using a vector of x coordinate points.
"""
function CylindricalGrid1D(
  x::AbstractVector{T},
  discretization_scheme::Symbol,
  snap_to_axis::Bool;
  backend=CPU(),
  on_bc=nothing,
  is_static=false,
  kwargs...,
) where {T}

  #

  scheme_name = Symbol(uppercase("$discretization_scheme"))
  nhalo = nhalo_lookup[scheme_name]

  if scheme_name === :MEG6 ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
  elseif scheme_name === :MEG6_SYMMETRIC ||
    discretization_scheme == :MontoneExplicitGradientScheme6thOrder
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder
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
    use_cache=true,
    celldims=size(domain_iterators.cell.full),
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=false,
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
    discr_scheme,
    _on_bc,
    snap_to_axis,
    is_static,
    scheme_name,
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
