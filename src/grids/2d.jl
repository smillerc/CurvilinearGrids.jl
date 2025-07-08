
"""
A 2D curvilinear grid. Curvilinear implies that the coorindates must be described
by 2 indices, e.g. `x -> x(i,j)` and `y -> y(i,j)`
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
  is_static::Bool
  is_orthogonal::Bool
  discretization_scheme_name::Symbol
end

"""
A 2D rectilinear grid. Rectilinear implies that the coorindates can be described
by 1D vectors, e.g. `x -> x(i)` and `y -> y(j)`
"""
struct RectilinearGrid2D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid2D
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
  is_static::Bool
  discretization_scheme_name::Symbol
end

"""
A 2D uniform grid. Grid coordinates can change over time, but transformations must be rigid
in order to maintain coordinate uniformity
"""
struct UniformGrid2D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid2D
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
  is_static::Bool
  discretization_scheme_name::Symbol
end

"""
AxisymmetricGrid2D
"""
struct AxisymmetricGrid2D{CO,CE,NV,EM,CM,DL,CI,DS} <: AbstractCurvilinearGrid2D
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
  snap_to_axis::Bool
  rotational_axis::Symbol
  is_static::Bool
  is_orthogonal::Bool
  discretization_scheme_name::Symbol
end

"""
    CurvilinearGrid2D((x0, y0), (x1, y1), (ni_cells, nj_cells)::NTuple{2,Int}, discretization_scheme::Symbol; T=Float64)

Create a `CurvilinearGrid2D` with start/end points `(x0, y0), (x1, y1)`, and the number of cells in each dimension. This will create a rectilinear-style
grid, but still remain curvilinear, such that the grid can dynamically change. For a grid that stays purely rectilinear, use the `RectilinearGrid2D` constructor instead.
"""
function CurvilinearGrid2D(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  discretization_scheme::Symbol;
  T=Float64,
  kwargs...,
)
  x = range(x0, x1; length=ni_cells + 1) .|> T
  y = range(y0, y1; length=nj_cells + 1) .|> T

  return CurvilinearGrid2D(x, y, discretization_scheme; T=T, kwargs...)
end

"""
    CurvilinearGrid2D(x::AbstractVector{T}, y::AbstractVector{T}, discretization_scheme::Symbol; kwargs...) where {T}

Create a `CurvilinearGrid2D` with 1D `x` and `y` coordinate vectors. The input coordinates do not include halo / ghost data since the geometry is undefined in these regions.
"""
function CurvilinearGrid2D(
  x::AbstractVector{T}, y::AbstractVector{T}, discretization_scheme::Symbol; kwargs...
) where {T}
  ni = length(x)
  nj = length(y)

  X = zeros(T, ni, nj)
  Y = zeros(T, ni, nj)
  for j in eachindex(y)
    for i in eachindex(x)
      X[i, j] = x[i]
      Y[i, j] = y[j]
    end
  end

  return CurvilinearGrid2D(X, Y, discretization_scheme; kwargs...)
end

"""
    CurvilinearGrid2D(x::AbstractArray{T,2}, y::AbstractArray{T,2}, discretization_scheme=::Symbol; backend=CPU()) where {T}

Create a `CurvilinearGrid2D` with 2D `x` and `y` coordinates. The input coordinates do not include halo / ghost data since the geometry is undefined in these regions.
"""
function CurvilinearGrid2D(
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=false,
  is_orthogonal=false,
  init_metrics=true,
  empty_metrics=false,
  kwargs...,
) where {T}
  m = CurvilinearGrid2D(
    _grid_constructor(
      x,
      y,
      :curvilinear,
      discretization_scheme;
      backend=backend,
      is_static=is_static,
      is_orthogonal=is_orthogonal,
      empty_metrics=empty_metrics,
    )...,
  )

  if init_metrics && !empty_metrics
    update!(m; force=true)
  end

  return m
end

"""
    RectilinearGrid2D(x::AbstractVector{T}, y::AbstractVector{T}, discretization_scheme=::Symbol; backend=CPU()) where {T}

Create a `RectilinearGrid2D` with `x` and `y` coordinate vectors. The input coordinates do not include halo / ghost data since the geometry is undefined in these regions.
The `RectilinearGrid2D` type uses specialized storage based on the assertion that x and y are defined purely by 1D vectors, rather than 2D arrays. Any dynamic motion must
obey this restriction. For a fully dynamic mesh with 2D varying coordinates, use a `CurvilinearGrid2D` type instead.
"""
function RectilinearGrid2D(
  x::AbstractVector{T},
  y::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=true,
  init_metrics=true,
  empty_metrics=false,
) where {T}
  ni = length(x)
  nj = length(y)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  if nj < 2
    error("The y vector must have more than 2 points")
  end

  if !all(diff(x) .> 0)
    error("Invalid x vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(y) .> 0)
    error("Invalid y vector, spacing between vertices must be > 0 everywhere")
  end

  x2d = zeros(T, ni, nj)
  y2d = zeros(T, ni, nj)

  @inbounds for j in 1:nj
    for i in 1:ni
      x2d[i, j] = x[i]
      y2d[i, j] = y[j]
    end
  end

  m = RectilinearGrid2D(
    _grid_constructor(
      x2d,
      y2d,
      :rectilinear,
      discretization_scheme;
      backend=backend,
      is_static=is_static,
      empty_metrics=empty_metrics,
    )...,
  )

  if init_metrics && !empty_metrics
    update!(m; force=true)
  end
  return m
end

"""
    RectilinearGrid2D((x0, y0), (x1, y1), (ni_cells, nj_cells)::NTuple{2,Int}, discretization_scheme::Symbol; backend=CPU(), is_static=true, init_metrics=true, empty_metrics=false)

Create a `RectilinearGrid2D` with start/end points `(x0, y0), (x1, y1)`, and the number of cells in each dimension.
The `RectilinearGrid2D` type uses specialized storage based on the assertion that x and y are defined purely by 1D vectors, rather than 2D arrays. Any dynamic motion must
obey this restriction. For a fully dynamic mesh with 2D varying coordinates, use a `CurvilinearGrid2D` type instead.
"""
function RectilinearGrid2D(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=true,
  init_metrics=true,
  empty_metrics=false,
)
  if ni_cells < 2 || nj_cells < 2
    error(
      "The number of cells specified must be > 1, given cell dims are $((ni_cells, nj_cells))",
    )
  end

  return RectilinearGrid2D(
    collect(range(x0, x1; length=ni_cells + 1)),
    collect(range(y0, y1; length=nj_cells + 1)),
    discretization_scheme;
    backend=backend,
    is_static=is_static,
    init_metrics=init_metrics,
    empty_metrics=empty_metrics,
  )
end

"""
    UniformGrid2D((x0, y0), (x1, y1), ∂x, discretization_scheme=::Symbol; backend=CPU()) where {T}

Create a 2d uniform grid with a grid spacing and a shape. The input coordinates do not include halo / ghost data since the geometry is undefined in these regions.
The `UniformGrid2D` type uses specialized storage based on the assertion that x and y uniform and can be stored with reduced date. Any dynamic motion must
obey this restriction. For a fully dynamic mesh with 2D varying coordinates, use a `CurvilinearGrid2D` type instead.
"""
function UniformGrid2D(
  (x0, y0),
  (x1, y1),
  ∂x::T,
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=true,
  init_metrics=true,
  empty_metrics=false,
) where {T<:Real}
  x = Vector(x0:∂x:x1)
  y = Vector(y0:∂x:y1)

  ni = length(x)
  nj = length(y)

  if ni < 2
    error("The x vector must have more than 2 points")
  end

  if nj < 2
    error("The y vector must have more than 2 points")
  end

  if !all(diff(x) .> 0)
    error("Invalid x vector, spacing between vertices must be > 0 everywhere")
  end

  if !all(diff(y) .> 0)
    error("Invalid y vector, spacing between vertices must be > 0 everywhere")
  end

  x2d = zeros(T, ni, nj)
  y2d = zeros(T, ni, nj)

  @inbounds for j in 1:nj
    for i in 1:ni
      x2d[i, j] = x[i]
      y2d[i, j] = y[j]
    end
  end

  m = UniformGrid2D(
    _grid_constructor(
      x2d,
      y2d,
      :uniform,
      discretization_scheme;
      backend=backend,
      is_static=is_static,
      empty_metrics=empty_metrics,
    )...,
  )

  if init_metrics && !empty_metrics
    update!(m; force=true)
  end
  return m
end

function _grid_constructor(
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  tag::Symbol,
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=false,
  is_orthogonal=false,
  empty_metrics=false,
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

  if tag == :rectilinear
    if empty_metrics
      cell_center_metrics, edge_metrics = (nothing, nothing)
    else
      cell_center_metrics, edge_metrics = get_metric_soa_rectilinear2d(celldims, backend, T)
    end
  elseif tag == :uniform
    if empty_metrics
      cell_center_metrics, edge_metrics = (nothing, nothing)
    else
      cell_center_metrics, edge_metrics = get_metric_soa_uniform2d(celldims, backend, T)
    end
  elseif tag == :curvilinear
    if empty_metrics
      cell_center_metrics, edge_metrics = (nothing, nothing)
    else
      cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)
    end
  end

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

  discr_scheme = MetricDiscretizationScheme(;
    use_cache=true,
    celldims=size(domain_iterators.cell.full),
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=false,
  )

  if tag == :rectilinear || tag === :uniform
    return (
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
      is_static,
      scheme_name,
    )
  elseif tag === :curvilinear
    return (
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
      is_static,
      is_orthogonal,
      scheme_name,
    )
  end
end

"""
    AxisymmetricGrid2D((x0, y0), (x1, y1), (ni_cells, nj_cells)::NTuple{2,Int}, discretization_scheme::Symbol, snap_to_axis::Bool, rotational_axis::Symbol; T=Float64)

Create an axisymmetric 2d grid with end points as `(x0, y0)`, and `(x1, y1)`, and `(ni_cells, nj_cells)`, with the symmetry axis `rotational_axis = :x` or `rotational_axis = :y`. 
Enforce coordinates to stay on the axis via `snap_to_axis=true`. The input coordinates do not include halo / ghost data since the geometry is undefined in these regions. 
"""
function AxisymmetricGrid2D(
  (x0, y0),
  (x1, y1),
  (ni_cells, nj_cells)::NTuple{2,Int},
  discretization_scheme::Symbol,
  snap_to_axis::Bool,
  rotational_axis::Symbol;
  T=Float64,
  kwargs...,
)
  x = range(x0, x1; length=ni_cells + 1)
  y = range(y0, y1; length=nj_cells + 1)

  ni = length(x)
  nj = length(y)

  X = zeros(T, ni, nj)
  Y = zeros(T, ni, nj)
  for j in eachindex(y)
    for i in eachindex(x)
      X[i, j] = x[i]
      Y[i, j] = y[j]
    end
  end

  return AxisymmetricGrid2D(
    X, Y, discretization_scheme, snap_to_axis, rotational_axis; kwargs...
  )
end

"""
    AxisymmetricGrid2D(x::AbstractMatrix{T}, y::AbstractMatrix{T},  nhalo::Int,  snap_to_axis::Bool,  rotational_axis::Symbol;  T=Float64, backend=CPU())

Create an axisymmetric 2d grid with `x` and `y` coordinates, with the symmetry axis `rotational_axis = :x` or `rotational_axis = :y`. Enforce coordinates to stay on the axis via `snap_to_axis=true`.
The input coordinates do not include halo / ghost data since the geometry is undefined in these regions. 
"""
function AxisymmetricGrid2D(
  x::AbstractMatrix{T},
  y::AbstractMatrix{T},
  discretization_scheme,
  snap_to_axis::Bool,
  rotational_axis::Symbol;
  backend=CPU(),
  is_static=false,
  is_orthogonal=false,
  make_uniform=false,
  kwargs...,
) where {T}

  #

  scheme_name = Symbol(uppercase("$discretization_scheme"))
  nhalo = nhalo_lookup[scheme_name]

  if (
    scheme_name === :MEG6 ||
    scheme_name === :MEG6_SYMMETRIC ||
    discretization_scheme === :MontoneExplicitGradientScheme6thOrder
  )
    MetricDiscretizationScheme = MontoneExplicitGradientScheme6thOrder

  else
    error("Only MontoneExplicitGradientScheme6thOrder or MEG6 is supported for now")
  end

  @assert size(x) == size(y)
  @assert rotational_axis === :x || rotational_axis === :y

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

  # _centroid_coordinates!(centroids, coords, domain_iterators.cell.domain)

  node_velocities = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  discr_scheme = MetricDiscretizationScheme(;
    use_cache=true,
    celldims=size(domain_iterators.cell.full),
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=false,
  )

  m = AxisymmetricGrid2D(
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
    snap_to_axis,
    rotational_axis,
    is_static,
    is_orthogonal,
    scheme_name,
  )

  update!(m)

  if make_uniform
    first_idx = first(m.iterators.cell.domain)

    fill!(m.cell_center_metrics.J, m.cell_center_metrics.J[first_idx])

    fill!(m.cell_center_metrics.x₁.ξ, m.cell_center_metrics.x₁.ξ[first_idx])
    fill!(m.cell_center_metrics.x₂.ξ, m.cell_center_metrics.x₂.ξ[first_idx])
    fill!(m.cell_center_metrics.x₁.η, m.cell_center_metrics.x₁.η[first_idx])
    fill!(m.cell_center_metrics.x₂.η, m.cell_center_metrics.x₂.η[first_idx])

    fill!(m.cell_center_metrics.ξ.x₁, m.cell_center_metrics.ξ.x₁[first_idx])
    fill!(m.cell_center_metrics.ξ.x₂, m.cell_center_metrics.ξ.x₂[first_idx])
    fill!(m.cell_center_metrics.η.x₁, m.cell_center_metrics.η.x₁[first_idx])
    fill!(m.cell_center_metrics.η.x₂, m.cell_center_metrics.η.x₂[first_idx])

    fill!(m.edge_metrics.i₊½.J, m.edge_metrics.i₊½.J[first_idx])

    fill!(m.edge_metrics.i₊½.ξ̂.x₁, m.edge_metrics.i₊½.ξ̂.x₁[first_idx])
    fill!(m.edge_metrics.i₊½.ξ̂.x₂, m.edge_metrics.i₊½.ξ̂.x₂[first_idx])
    fill!(m.edge_metrics.i₊½.ξ̂.t, m.edge_metrics.i₊½.ξ̂.t[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.x₁, m.edge_metrics.i₊½.η̂.x₁[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.x₂, m.edge_metrics.i₊½.η̂.x₂[first_idx])
    fill!(m.edge_metrics.i₊½.η̂.t, m.edge_metrics.i₊½.η̂.t[first_idx])

    fill!(m.edge_metrics.j₊½.J, m.edge_metrics.j₊½.J[first_idx])

    fill!(m.edge_metrics.j₊½.ξ̂.x₁, m.edge_metrics.j₊½.ξ̂.x₁[first_idx])
    fill!(m.edge_metrics.j₊½.ξ̂.x₂, m.edge_metrics.j₊½.ξ̂.x₂[first_idx])
    fill!(m.edge_metrics.j₊½.ξ̂.t, m.edge_metrics.j₊½.ξ̂.t[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.x₁, m.edge_metrics.j₊½.η̂.x₁[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.x₂, m.edge_metrics.j₊½.η̂.x₂[first_idx])
    fill!(m.edge_metrics.j₊½.η̂.t, m.edge_metrics.j₊½.η̂.t[first_idx])
  end

  return m
end

"""Update metrics after grid coordinates change"""
function update!(mesh::AbstractCurvilinearGrid2D; force=false)
  if !mesh.is_static || force
    _centroid_coordinates!(
      mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
    )
    update_metrics!(mesh)
    _check_valid_metrics(mesh)
  else
    error("Attempting to update grid metrics when grid.is_static = true!")
  end
  return nothing
end

"""Update metrics after grid coordinates change"""
function update!(mesh::AxisymmetricGrid2D)
  _centroid_coordinates!(
    mesh.centroid_coordinates, mesh.node_coordinates, mesh.iterators.cell.domain
  )
  update_metrics!(mesh)

  if mesh.snap_to_axis
    _snap_nodes_to_axis(mesh)
    # else
    #   _check_nodes_along_axis(mesh)
  end

  _check_valid_metrics(mesh)
  return nothing
end

function jacobian_matrix(mesh::AbstractCurvilinearGrid2D, (i, j))
  xξ = mesh.cell_center_metrics.x₁.ξ
  yξ = mesh.cell_center_metrics.x₂.ξ
  xη = mesh.cell_center_metrics.x₁.η
  yη = mesh.cell_center_metrics.x₂.η

  return @SMatrix [
    xξ[i, j] xη[i, j]
    yξ[i, j] yη[i, j]
  ]
end

"""
    jacobian(mesh::CurvilinearGrid2D, idx)
    jacobian(mesh::RectilinearGrid2D, idx)
    jacobian(mesh::UniformGrid2D, idx)

The cell-centroid Jacobian (determinant of the Jacobian matrix)
"""
jacobian(mesh::CurvilinearGrid2D, idx) = det(jacobian_matrix(mesh, idx))

jacobian(mesh::RectilinearGrid2D, idx) = det(jacobian_matrix(mesh, idx))

jacobian(mesh::UniformGrid2D, idx) = det(jacobian_matrix(mesh, idx))

"""
    forward_cell_metrics(mesh::AbstractCurvilinearGrid2D, idx)

Get the forward cell metrics for a given cell index `idx`
"""
function forward_cell_metrics(mesh::AbstractCurvilinearGrid2D, idx)
  (; x=mesh.cell_center_metrics.x₁[idx], y=mesh.cell_center_metrics.x₂[idx])
end

"""
    inverse_cell_metrics(mesh::AbstractCurvilinearGrid2D, idx)

Get the inverse cell metrics for a given cell index `idx`
"""
function inverse_cell_metrics(mesh::AbstractCurvilinearGrid2D, idx)
  (;
    ξ=mesh.cell_center_metrics.ξ[idx],
    η=mesh.cell_center_metrics.η[idx],
    ξ̂=mesh.cell_center_metrics.ξ̂[idx],
    η̂=mesh.cell_center_metrics.η̂[idx],
  )
end

# ------------------------------------------------------------------
# Velocity Functions
# ------------------------------------------------------------------

@inline grid_velocities(::AbstractCurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0) =
  (0.0, 0.0)
# @inline centroid_velocities(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)
# @inline node_velocities(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t) = (0.0, 0.0)

# ------------------------------------------------------------------
# Private Functions
# ------------------------------------------------------------------

function _centroid_coordinates!(
  centroids::StructArray{T,2}, coords::StructArray{T,2}, domain
) where {T}
  x = coords.x
  y = coords.y
  # Populate the centroid coordinates
  for idx in domain
    i, j = idx.I
    centroids.x[idx] = 0.25(x[i, j] + x[i + 1, j] + x[i + 1, j + 1] + x[i, j + 1])
    centroids.y[idx] = 0.25(y[i, j] + y[i + 1, j] + y[i + 1, j + 1] + y[i, j + 1])
  end

  return nothing
end

function _check_valid_metrics(mesh::AbstractCurvilinearGrid2D)
  domain = mesh.iterators.cell.domain
  i₊½_domain = expand(domain, 1, -1)
  j₊½_domain = expand(domain, 2, -1)

  @views begin
    centroid_metrics_valid =
      all(isfinite.(mesh.cell_center_metrics.J[domain])) &&
      all(mesh.cell_center_metrics.J[domain] .> 0) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.ξ.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₁[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.η.x₂[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₁.η[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.ξ[domain])) &&
      all(isfinite.(mesh.cell_center_metrics.x₂.η[domain]))

    edge_metrics_valid =
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.ξ̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₁[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.x₂[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.i₊½.η̂.t[i₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.ξ̂.t[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₁[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.x₂[j₊½_domain])) &&
      all(isfinite.(mesh.edge_metrics.j₊½.η̂.t[j₊½_domain]))
  end

  if !edge_metrics_valid
    error("Invalid edge metrics found")
  end

  if !centroid_metrics_valid
    error("Invalid centroid metrics found")
  end

  return nothing
end

function _check_nodes_along_axis(mesh::AxisymmetricGrid2D)
  domain = mesh.iterators.node.domain

  @views begin
    if mesh.rotational_axis === :x
      axis_domain = domain[:, 1]

      if any(!iszero(mesh.node_coordinates.y[axis_domain]))
        error(
          "Nodes not aligned to axis of symmetry (`snap_to_axis = true`). Set snap_to_axis = false to disable checks",
        )
      end
    else # mesh.rotational_axis === :y
      axis_domain = domain[1, :]

      if any(!iszero(mesh.node_coordinates.x[axis_domain]))
        error(
          "Nodes not aligned to axis of symmetry (`snap_to_axis = true`). Set snap_to_axis = false to disable checks",
        )
      end
    end
  end
end

function _snap_nodes_to_axis(mesh::AxisymmetricGrid2D)
  domain = mesh.iterators.node.domain

  @views begin
    if mesh.rotational_axis === :x
      axis_domain = domain[:, 1]
      mesh.node_coordinates.y[axis_domain] .= 0

    else # mesh.rotational_axis === :y
      axis_domain = domain[1, :]
      mesh.node_coordinates.x[axis_domain] .= 0
    end
  end
end
