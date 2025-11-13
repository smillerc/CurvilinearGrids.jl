
struct CurvilinearGrid3D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid3D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{3,Int}
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  is_static::Bool
  is_orthogonal::Bool
  discretization_scheme_name::Symbol
  halo_coords_included::Bool
end

struct RectilinearGrid3D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid3D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{3,Int}
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  is_static::Bool
  discretization_scheme_name::Symbol
  halo_coords_included::Bool
end

struct UniformGrid3D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid3D
  node_coordinates::CO
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{3,Int}
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  is_static::Bool
  discretization_scheme_name::Symbol
  halo_coords_included::Bool
end

"""
    CurvilinearGrid3D(x, y, z, discretization_scheme::Symbol; backend=CPU(), is_static=false, is_orthogonal=false, tiles=nothing, init_metrics=true)

Construct a curvilinear grid in 3D using 3D arrays of x/y/z coordinates.
"""
function CurvilinearGrid3D(
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=false,
  is_orthogonal=false,
  init_metrics=true,
  empty_metrics=false,
  halo_coords_included=false,
  kwargs...,
) where {T}

  #
  GradientDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )

  # limits, iterators = get_node_cell_iterators(x, y, z, nhalo)
  limits, iterators = get_iterators(size(x), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = GradientDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes = _grid_constructor(
    x, y, z, iterators, halo_coords_included; backend=backend
  )
  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  m = CurvilinearGrid3D(
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
    is_orthogonal,
    scheme_name,
    halo_coords_included,
  )

  if init_metrics && !empty_metrics
    update!(m, true, halo_coords_included)
  end

  return m
end

"""
    CurvilinearGrid3D(x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, discretization_scheme::Symbol; kwargs...)

Construct a curvilinear grid in 3D using x/y/z vector coordinates.
"""
function CurvilinearGrid3D(
  x::AbstractVector{T},
  y::AbstractVector{T},
  z::AbstractVector{T},
  discretization_scheme::Symbol;
  kwargs...,
) where {T}
  ni = length(x)
  nj = length(y)
  nk = length(z)

  x3d = zeros(T, ni, nj, nk)
  y3d = zeros(T, ni, nj, nk)
  z3d = zeros(T, ni, nj, nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x3d[i, j, k] = x[i]
        y3d[i, j, k] = y[j]
        z3d[i, j, k] = z[k]
      end
    end
  end

  return CurvilinearGrid3D(x3d, y3d, z3d, discretization_scheme; kwargs...)
end

"""
    CurvilinearGrid3D((x0, y0, z0), (x1, y1, z1), (ni_cells, nj_cells, nk_cells)::NTuple{3,Int}, discretization_scheme::Symbol; kwargs...)

Construct a curvilinear grid in 3D using start/end coordinates and dimensions along each axis.
"""
function CurvilinearGrid3D(
  (x0, y0, z0),
  (x1, y1, z1),
  (ni_cells, nj_cells, nk_cells)::NTuple{3,Int},
  discretization_scheme::Symbol;
  T=Float64,
  kwargs...,
)
  return CurvilinearGrid3D(
    range(x0, x1; length=ni_cells + 1) .|> T,
    range(y0, y1; length=nj_cells + 1) .|> T,
    range(z0, z1; length=nk_cells + 1) .|> T,
    discretization_scheme;
    kwargs...,
  )
end

"""
    RectilinearGrid3D(x, y, z, discretization_scheme::Symbol; backend=CPU(), is_static=false, is_static=true,init_metrics=true)
    RectilinearGrid3D((x0, y0, z0), (x1, y1, z1), (∂x, ∂y, ∂z)::NTuple{3,AbstractFloat}, discretization_scheme::Symbol; backend=CPU(), is_static=true, init_metrics=true)
    RectilinearGrid3D((x0, y0, z0), (x1, y1, z1), (ni, nj, nk)::NTuple{3,Int}, discretization_scheme::Symbol; backend=CPU(), is_static=true, init_metrics=true)

Construct a rectilinear grid in 3D using 3D arrays of x/y/z coordinates. This constructor utilizes `RectilinearArray`s to optimize 
data storage, and therefore should only be used with grids that are rectilinear.
"""
function RectilinearGrid3D(
  x::AbstractVector{T},
  y::AbstractVector{T},
  z::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=true,
  init_metrics=true,
  empty_metrics=false,
  halo_coords_included=false,
  kwargs...,
) where {T}

  #
  ni = length(x)
  nj = length(y)
  nk = length(z)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  if nj < 2
    error("The y vector must have more than 2 points")
  end

  if nk < 2
    error("The z vector must have more than 2 points")
  end

  if !all(diff(x) .> 0)
    error("Invalid x vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(y) .> 0)
    error("Invalid y vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(z) .> 0)
    error("Invalid z vector, spacing between vertices must be > 0 everywhere")
  end

  x3d = zeros(T, ni, nj, nk)
  y3d = zeros(T, ni, nj, nk)
  z3d = zeros(T, ni, nj, nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x3d[i, j, k] = x[i]
        y3d[i, j, k] = y[j]
        z3d[i, j, k] = z[k]
      end
    end
  end

  #
  GradientDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )

  # limits, iterators = get_node_cell_iterators(x3d, y3d, z3d, nhalo)
  limits, iterators = get_iterators(size(x3d), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = GradientDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes = _grid_constructor(
    x3d, y3d, z3d, iterators, halo_coords_included; backend=backend
  )

  cell_center_metrics, edge_metrics = get_metric_soa_rectilinear3d(celldims, backend, T)

  m = RectilinearGrid3D(
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
    halo_coords_included,
  )

  if init_metrics && !empty_metrics
    update!(m, true, halo_coords_included)
  end

  return m
end

function RectilinearGrid3D(
  (x0, y0, z0),
  (x1, y1, z1),
  (ni_cells, nj_cells, nk_cells)::NTuple{3,Int},
  discretization_scheme::Symbol;
  kwargs...,
)
  return RectilinearGrid3D(
    range(x0, x1; length=ni_cells + 1),
    range(y0, y1; length=nj_cells + 1),
    range(z0, z1; length=nk_cells + 1),
    discretization_scheme;
    kwargs...,
  )
end

"""
    UniformGrid3D((x0, y0, z0), (x1, y1, z1), ∂x, shape::NTuple{3, Int}, discretization_scheme::Symbol; backend=CPU(), is_static=false,is_static=true, init_metrics=true)

Construct a uniform grid in 3D using a grid spacing and a shape tuple. This constructor utilizes `RectilinearArray`s to 
optimize data storage, and therefore should only be used with grids that are uniform."""
function UniformGrid3D(
  (x0, y0, z0),
  (x1, y1, z1),
  ∂x::T,
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=true,
  init_metrics=true,
  empty_metrics=false,
  halo_coords_included=false,
  kwargs...,
) where {T<:Real}

  #
  x = x0:∂x:x1
  y = y0:∂x:y1
  z = z0:∂x:z1

  ni = length(x)
  nj = length(y)
  nk = length(z)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  if nj < 2
    error("The y vector must have more than 2 points")
  end

  if nk < 2
    error("The z vector must have more than 2 points")
  end

  if !all(diff(x) .> 0)
    error("Invalid x vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(y) .> 0)
    error("Invalid y vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(z) .> 0)
    error("Invalid z vector, spacing between vertices must be > 0 everywhere")
  end

  x3d = zeros(T, ni, nj, nk)
  y3d = zeros(T, ni, nj, nk)
  z3d = zeros(T, ni, nj, nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        x3d[i, j, k] = x[i]
        y3d[i, j, k] = y[j]
        z3d[i, j, k] = z[k]
      end
    end
  end
  #
  GradientDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )

  # limits, iterators = get_node_cell_iterators(x3d, y3d, z3d, nhalo)
  limits, iterators = get_iterators(size(x3d), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = GradientDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes = _grid_constructor(
    x3d, y3d, z3d, iterators, halo_coords_included; backend=backend
  )

  cell_center_metrics, edge_metrics = get_metric_soa_uniform3d(celldims, backend, T)

  m = UniformGrid3D(
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
    halo_coords_included,
  )

  if init_metrics && !empty_metrics
    update!(m, true, halo_coords_included)
  end

  return m
end

# common grid constructor for uniform grids
function _grid_constructor(
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  domain_iterators,
  halo_coords_included;
  backend=KernelAbstractions.CPU(),
  kwargs...,
) where {T}

  #
  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  centroids = StructArray((
    x=KernelAbstractions.zeros(backend, T, celldims),
    y=KernelAbstractions.zeros(backend, T, celldims),
    z=KernelAbstractions.zeros(backend, T, celldims),
  ))

  if halo_coords_included
    copy!(coords.x, x)
    copy!(coords.y, y)
    copy!(coords.z, z)
  else
    @views begin
      copy!(coords.x[domain_iterators.node.domain], x)
      copy!(coords.y[domain_iterators.node.domain], y)
      copy!(coords.z[domain_iterators.node.domain], z)
    end
  end

  if halo_coords_included
    domain = domain_iterators.cell.full
  else
    domain = domain_iterators.cell.domain
  end

  _centroid_coordinates_kernel!(backend)(centroids, coords, domain; ndrange=size(domain))

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  nnodes = size(x)
  return (coords, centroids, node_velocities, nnodes)
end

"""Update metrics after grid coordinates change"""
function update!(mesh::AbstractCurvilinearGrid3D, force::Bool, include_halo_region::Bool)
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
      metric_domain;
      ndrange=size(metric_domain),
    )
    update_metrics!(mesh; include_halo_region=include_halo_region)
    _check_valid_metrics(mesh)
  else
    @warn("Attempting to update grid metrics when grid.is_static = true!")
  end
  return nothing
end

function jacobian_matrix(mesh::CurvilinearGrid3D, (i, j, k))
  xξ = mesh.cell_center_metrics.x₁.ξ
  yξ = mesh.cell_center_metrics.x₂.ξ
  zξ = mesh.cell_center_metrics.x₂.ξ
  xη = mesh.cell_center_metrics.x₁.η
  yη = mesh.cell_center_metrics.x₂.η
  zη = mesh.cell_center_metrics.x₂.η
  xζ = mesh.cell_center_metrics.x₁.η
  yζ = mesh.cell_center_metrics.x₂.η
  zζ = mesh.cell_center_metrics.x₂.η

  return @SMatrix [
    xξ[i, j, k] xη[i, j, k] xζ[i, j, k]
    yξ[i, j, k] yη[i, j, k] yζ[i, j, k]
    zξ[i, j, k] zη[i, j, k] zζ[i, j, k]
  ]
end

"""
    forward_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)

Get the forward cell metrics for a given cell index `idx`
"""
function forward_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)
  (;
    x=mesh.cell_center_metrics.x₁[idx],
    y=mesh.cell_center_metrics.x₂[idx],
    z=mesh.cell_center_metrics.x₃[idx],
  )
end

"""
    inverse_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)

Get the inverse cell metrics for a given cell index `idx`
"""
function inverse_cell_metrics(mesh::AbstractCurvilinearGrid3D, idx)
  (;
    ξ=mesh.cell_center_metrics.ξ[idx],
    η=mesh.cell_center_metrics.η[idx],
    ζ=mesh.cell_center_metrics.ζ[idx],
    ξ̂=mesh.cell_center_metrics.ξ̂[idx],
    η̂=mesh.cell_center_metrics.η̂[idx],
    ζ̂=mesh.cell_center_metrics.ζ̂[idx],
  )
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0) =
  (0.0, 0.0, 0.0)
# @inline centroid_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t) = (0.0, 0.0, 0.0)
# @inline node_velocities(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t) = (0.0, 0.0, 0.0)

# ------------------------------------------------------------------
# Private Functions
# ------------------------------------------------------------------

@kernel inbounds = true function _centroid_coordinates_kernel!(
  centroids::StructArray{T,3}, coords::StructArray{T,3}, domain
) where {T}
  idx = @index(Global, Linear)
  didx = domain[idx]
  i, j, k = didx.I

  x = coords.x
  y = coords.y
  z = coords.z

  centroids.x[didx] =
    0.125(
      x[i, j, k] +
      x[i + 1, j, k] +
      x[i + 1, j + 1, k] +
      x[i, j + 1, k] +
      x[i, j, k + 1] +
      x[i + 1, j, k + 1] +
      x[i + 1, j + 1, k + 1] +
      x[i, j + 1, k + 1]
    )

  centroids.y[didx] =
    0.125(
      y[i, j, k] +
      y[i + 1, j, k] +
      y[i + 1, j + 1, k] +
      y[i, j + 1, k] +
      y[i, j, k + 1] +
      y[i + 1, j, k + 1] +
      y[i + 1, j + 1, k + 1] +
      y[i, j + 1, k + 1]
    )

  centroids.z[didx] =
    0.125(
      z[i, j, k] +
      z[i + 1, j, k] +
      z[i + 1, j + 1, k] +
      z[i, j + 1, k] +
      z[i, j, k + 1] +
      z[i + 1, j, k + 1] +
      z[i + 1, j + 1, k + 1] +
      z[i, j + 1, k + 1]
    )
end

function _check_valid_metrics(mesh::AbstractCurvilinearGrid3D)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, 1, -1)
  j₊½_domain = expand(domain, 2, -1)
  k₊½_domain = expand(domain, 3, -1)

  @views begin
    centroid_metrics_valid =
      all(isfinite.(mesh.cell_center_metrics.J[domain])) &&
      all(mesh.cell_center_metrics.J[domain] .> 0) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₃[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₃[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ζ.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ζ.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ζ.x₃[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.η[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.ζ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.η[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.ζ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₃.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₃.η[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₃.ζ[domain]))

    edge_metrics_valid =
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₃[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₃[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.x₃[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ζ̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₃[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.t[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₃[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.t[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.x₃[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ζ̂.t[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.x₁[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.x₂[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.x₃[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ξ̂.t[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.x₁[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.x₂[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.x₃[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.η̂.t[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.x₁[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.x₂[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.x₃[k₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.k₊½.ζ̂.t[k₊½_domain]))
  end

  if !edge_metrics_valid
    error("Invalid edge metrics found")
  end

  if !centroid_metrics_valid
    error("Invalid centroid metrics found")
  end

  return nothing
end
