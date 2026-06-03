"""
Unified grid model (Phase 1, additive):
- MappedGrid
- DiscreteGrid
- OrthogonalGrid
"""

"""
Base abstract type for all unified-grid implementations.
"""
abstract type AbstractUnifiedGrid end

"""
Base abstract type for unified grids that own mapping/metric caches.
"""
abstract type AbstractMappedOrDiscreteGrid <: AbstractUnifiedGrid end

#
# Independent metric caches
#

"""
    UnifiedMetricCache{D}

Cache wrapper for one metric domain (cell or face).

# Fields
  - `data`: Backing metric storage container.
  - `valid`: Cache validity flag.
  - `mode`: Cache mode (`:eager`, `:lazy`, `:off`).
"""
mutable struct UnifiedMetricCache{D}
  data::D
  valid::Bool
  mode::Symbol
end

"""
    UnifiedMetricCaches{C,F}

Container for independently managed cell and face metric caches.

# Fields
  - `cell`: Cell-center metric cache.
  - `face`: Face metric cache.
"""
mutable struct UnifiedMetricCaches{C,F}
  cell::UnifiedMetricCache{C}
  face::UnifiedMetricCache{F}
end

function _check_cache_mode(cache_mode::Symbol)
  if cache_mode ∉ (:eager, :lazy, :off)
    throw(
      ArgumentError(
        "Invalid cache mode `$cache_mode`. Expected one of `:eager`, `:lazy`, `:off`."
      ),
    )
  end
  return cache_mode
end

function _new_metric_caches(cache_mode::Symbol, cell_data, face_data)
  mode = _check_cache_mode(cache_mode)
  UnifiedMetricCaches(
    UnifiedMetricCache(cell_data, false, mode), UnifiedMetricCache(face_data, false, mode)
  )
end

@inline function _has_metric_storage(grid::AbstractMappedOrDiscreteGrid)
  return grid.metric_caches !== nothing
end

@inline function _require_metric_storage(
  grid::AbstractMappedOrDiscreteGrid, caller::AbstractString
)
  if !_has_metric_storage(grid)
    throw(
      ArgumentError(
        "`$caller` is unavailable because this grid was created without metric storage (`compute_metrics=false, cache_mode=:off`).",
      ),
    )
  end
  return nothing
end

@inline function _has_metric_functions(grid::AbstractMappedOrDiscreteGrid)
  return grid.metric_functions_cache !== nothing
end

@inline function _require_metric_functions(
  grid::AbstractMappedOrDiscreteGrid, caller::AbstractString
)
  if !_has_metric_functions(grid)
    throw(
      ArgumentError(
        "`$caller` is unavailable because this grid has no metric function cache."
      ),
    )
  end
  return nothing
end

#
# Unified-grid geometry/metric helpers (AoS)
#

function _allocate_unified_cell_metric_storage(
  ::Val{N}, backend, ::Type{T}, iterators
) where {N,T}
  celldims = size(iterators.cell.full)
  metric_type = _metric_eltype(Val(N), T)
  metric_array() = KernelAbstractions.zeros(backend, metric_type, celldims...; unified=false)
  return (; forward=metric_array(), inverse=metric_array())
end

function _allocate_unified_face_metric_storage(
  ::Val{N}, backend, ::Type{T}, iterators
) where {N,T}
  celldims = size(iterators.cell.full)
  metric_type = _metric_eltype(Val(N), T)
  conserved_metric_type = _conserved_metric_eltype(Val(N), T)
  metric_array() = KernelAbstractions.zeros(backend, metric_type, celldims...; unified=false)
  function conserved_metric_array()
    KernelAbstractions.zeros(backend, conserved_metric_type, celldims...; unified=false)
  end
  return ntuple(
    _ -> (;
      forward=metric_array(), inverse=metric_array(), conserved=conserved_metric_array()
    ),
    N,
  )
end

function _allocate_unified_metric_storage(
  ::Val{N}, backend, ::Type{T}, iterators
) where {N,T}
  cell = _allocate_unified_cell_metric_storage(Val(N), backend, T, iterators)
  face = _allocate_unified_face_metric_storage(Val(N), backend, T, iterators)

  return cell, face
end

function _ensure_metric_storage!(
  grid::AbstractMappedOrDiscreteGrid, ::Val{N}, ::Type{T}
) where {N,T}
  _require_metric_storage(grid, "_ensure_metric_storage!")
  return grid.metric_caches.cell.data, grid.metric_caches.face.data
end

@inline function _check_unified_basis_trait(basis::BasisTrait)
  if basis isa CartesianBasis || basis isa SphericalBasis
    return basis
  end
  throw(
    ArgumentError(
      "Unsupported basis trait $(typeof(basis)) for unified grids. Use `CartesianBasis()` or `SphericalBasis()`.",
    ),
  )
end

function _allocate_unified_coordinates(::Val{1}, iterators, backend, ::Type{T}) where {T}
  ni_nodes, = size(iterators.node.full)
  ni_cells, = size(iterators.cell.full)

  node_coordinates = (KernelAbstractions.zeros(backend, T, (ni_nodes,); unified=false),)
  centroid_coordinates = (KernelAbstractions.zeros(backend, T, (ni_cells,); unified=false),)
  face_coordinates = ((KernelAbstractions.zeros(backend, T, (ni_cells,); unified=false),),)
  return node_coordinates, centroid_coordinates, face_coordinates
end

function _allocate_unified_coordinates(::Val{2}, iterators, backend, ::Type{T}) where {T}
  ni_nodes, nj_nodes = size(iterators.node.full)
  ni_cells, nj_cells = size(iterators.cell.full)

  node_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes); unified=false),
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes); unified=false),
  )
  centroid_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells); unified=false),
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells); unified=false),
  )
  face_coordinates = ntuple(
    _ -> (
      KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells); unified=false),
      KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells); unified=false),
    ),
    Val(2),
  )
  return node_coordinates, centroid_coordinates, face_coordinates
end

function _allocate_unified_coordinates(::Val{3}, iterators, backend, ::Type{T}) where {T}
  ni_nodes, nj_nodes, nk_nodes = size(iterators.node.full)
  ni_cells, nj_cells, nk_cells = size(iterators.cell.full)

  node_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes); unified=false),
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes); unified=false),
    KernelAbstractions.zeros(backend, T, (ni_nodes, nj_nodes, nk_nodes); unified=false),
  )
  centroid_coordinates = (
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells); unified=false),
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells); unified=false),
    KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells); unified=false),
  )
  face_coordinates = ntuple(
    _ -> (
      KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells); unified=false),
      KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells); unified=false),
      KernelAbstractions.zeros(backend, T, (ni_cells, nj_cells, nk_cells); unified=false),
    ),
    Val(3),
  )
  return node_coordinates, centroid_coordinates, face_coordinates
end

@inline function _metric_cache_for_mapping(
  ::Val{1},
  mapping_functions,
  diff_backend,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  MetricCache(
    mapping_functions.x1, diff_backend; edge_interpolation_scheme=edge_interpolation_scheme
  )
end
@inline function _metric_cache_for_mapping(
  ::Val{2},
  mapping_functions,
  diff_backend,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  MetricCache(
    mapping_functions.x1,
    mapping_functions.x2,
    diff_backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )
end
@inline function _metric_cache_for_mapping(
  ::Val{3},
  mapping_functions,
  diff_backend,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
)
  MetricCache(
    mapping_functions.x1,
    mapping_functions.x2,
    mapping_functions.x3,
    diff_backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )
end

function _build_unified_components(
  ::Val{N},
  mapping_functions,
  celldims::NTuple{N,Int},
  nhalo::Int,
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait,
  backend,
  diff_backend,
  ::Type{T};
  global_cell_indices=nothing,
  build_metric_storage::Bool=true,
) where {N,T}
  iterators = get_iterators(celldims, nhalo, global_cell_indices)
  node_coordinates, centroid_coordinates, face_coordinates = _allocate_unified_coordinates(
    Val(N), iterators, backend, T
  )
  metric_functions_cache = _metric_cache_for_mapping(
    Val(N), mapping_functions, diff_backend, edge_interpolation_scheme
  )
  cell_metric_storage, face_metric_storage = if build_metric_storage
    _allocate_unified_metric_storage(Val(N), backend, T, iterators)
  else
    nothing, nothing
  end

  return (;
    node_coordinates,
    centroid_coordinates,
    face_coordinates,
    metric_functions_cache,
    cell_metric_storage,
    face_metric_storage,
    nhalo,
    iterators,
  )
end

@kernel function _compute_unified_nodes_1d_kernel!(
  x, x1_map, local_domain, global_domain, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  ξ = Iglobal.I[1] - nhalo
  x[I] = x1_map(t, ξ, params)
end

@kernel function _compute_unified_nodes_2d_kernel!(
  x, y, x1_map, x2_map, local_domain, global_domain, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  ξ = Iglobal.I[1] - nhalo
  η = Iglobal.I[2] - nhalo
  x[I] = x1_map(t, ξ, η, params)
  y[I] = x2_map(t, ξ, η, params)
end

@kernel function _compute_unified_nodes_3d_kernel!(
  x, y, z, x1_map, x2_map, x3_map, local_domain, global_domain, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  ξ = Iglobal.I[1] - nhalo
  η = Iglobal.I[2] - nhalo
  ζ = Iglobal.I[3] - nhalo
  x[I] = x1_map(t, ξ, η, ζ, params)
  y[I] = x2_map(t, ξ, η, ζ, params)
  z[I] = x3_map(t, ξ, η, ζ, params)
end

@kernel function _compute_unified_centroids_1d_kernel!(
  x, x1_map, local_domain, global_domain, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = one(eltype(x)) / 2
  ξ = Iglobal.I[1] - nhalo + half
  x[I] = x1_map(t, ξ, params)
end

@kernel function _compute_unified_centroids_2d_kernel!(
  x, y, x1_map, x2_map, local_domain, global_domain, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = one(eltype(x)) / 2
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  x[I] = x1_map(t, ξ, η, params)
  y[I] = x2_map(t, ξ, η, params)
end

@kernel function _compute_unified_centroids_3d_kernel!(
  x, y, z, x1_map, x2_map, x3_map, local_domain, global_domain, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = one(eltype(x)) / 2
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half
  x[I] = x1_map(t, ξ, η, ζ, params)
  y[I] = x2_map(t, ξ, η, ζ, params)
  z[I] = x3_map(t, ξ, η, ζ, params)
end

@kernel function _compute_unified_face_coordinates_1d_kernel!(
  x, x1_map, local_domain, global_domain, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  ξ = Iglobal.I[1] - nhalo + one(eltype(x))
  x[I] = x1_map(t, ξ, params)
end

@kernel function _compute_unified_face_coordinates_2d_kernel!(
  x, y, x1_map, x2_map, local_domain, global_domain, axis::Int, nhalo::Int, t, params
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = one(eltype(x)) / 2
  ξ = Iglobal.I[1] - nhalo + half + (axis == 1 ? half : zero(half))
  η = Iglobal.I[2] - nhalo + half + (axis == 2 ? half : zero(half))
  x[I] = x1_map(t, ξ, η, params)
  y[I] = x2_map(t, ξ, η, params)
end

@kernel function _compute_unified_face_coordinates_3d_kernel!(
  x,
  y,
  z,
  x1_map,
  x2_map,
  x3_map,
  local_domain,
  global_domain,
  axis::Int,
  nhalo::Int,
  t,
  params,
)
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = one(eltype(x)) / 2
  ξ = Iglobal.I[1] - nhalo + half + (axis == 1 ? half : zero(half))
  η = Iglobal.I[2] - nhalo + half + (axis == 2 ? half : zero(half))
  ζ = Iglobal.I[3] - nhalo + half + (axis == 3 ? half : zero(half))
  x[I] = x1_map(t, ξ, η, ζ, params)
  y[I] = x2_map(t, ξ, η, ζ, params)
  z[I] = x3_map(t, ξ, η, ζ, params)
end

function _compute_unified_node_coordinates!(
  node_coordinates, mapping_functions, iterators, nhalo::Int, t, params, ::Val{N}, backend
) where {N}
  _compute_unified_node_coordinates!(
    backend, node_coordinates, mapping_functions, iterators, nhalo, t, params, Val(N)
  )
end

function _compute_unified_node_coordinates!(
  backend::KernelAbstractions.CPU,
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  local_domain = iterators.node.full
  global_domain = iterators.global_domain.node.full
  _compute_unified_nodes_1d_kernel!(backend)(
    node_coordinates[1],
    mapping_functions.x1,
    local_domain,
    global_domain,
    nhalo,
    t,
    params;
    ndrange=size(local_domain),
  )
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_node_coordinates!(
  backend::KernelAbstractions.CPU,
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  local_domain = iterators.node.full
  global_domain = iterators.global_domain.node.full
  _compute_unified_nodes_2d_kernel!(backend)(
    node_coordinates[1],
    node_coordinates[2],
    mapping_functions.x1,
    mapping_functions.x2,
    local_domain,
    global_domain,
    nhalo,
    t,
    params;
    ndrange=size(local_domain),
  )
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_node_coordinates!(
  backend::KernelAbstractions.CPU,
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  local_domain = iterators.node.full
  global_domain = iterators.global_domain.node.full
  _compute_unified_nodes_3d_kernel!(backend)(
    node_coordinates[1],
    node_coordinates[2],
    node_coordinates[3],
    mapping_functions.x1,
    mapping_functions.x2,
    mapping_functions.x3,
    local_domain,
    global_domain,
    nhalo,
    t,
    params;
    ndrange=size(local_domain),
  )
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_node_coordinates!(
  backend::KernelAbstractions.Backend,
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  local_domain = iterators.node.full
  global_domain = iterators.global_domain.node.full
  T = eltype(node_coordinates[1])
  x_h = Array{T}(undef, size(node_coordinates[1]))

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξ = T(Iglobal.I[1] - nhalo)
    x_h[I] = mapping_functions.x1(t, ξ, params)
  end

  copyto!(node_coordinates[1], x_h)
  return nothing
end

function _compute_unified_node_coordinates!(
  backend::KernelAbstractions.Backend,
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  local_domain = iterators.node.full
  global_domain = iterators.global_domain.node.full
  T = eltype(node_coordinates[1])
  x_h = Array{T}(undef, size(node_coordinates[1]))
  y_h = Array{T}(undef, size(node_coordinates[2]))

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξ = T(Iglobal.I[1] - nhalo)
    η = T(Iglobal.I[2] - nhalo)
    x_h[I] = mapping_functions.x1(t, ξ, η, params)
    y_h[I] = mapping_functions.x2(t, ξ, η, params)
  end

  copyto!(node_coordinates[1], x_h)
  copyto!(node_coordinates[2], y_h)
  return nothing
end

function _compute_unified_node_coordinates!(
  backend::KernelAbstractions.Backend,
  node_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  local_domain = iterators.node.full
  global_domain = iterators.global_domain.node.full
  T = eltype(node_coordinates[1])
  x_h = Array{T}(undef, size(node_coordinates[1]))
  y_h = Array{T}(undef, size(node_coordinates[2]))
  z_h = Array{T}(undef, size(node_coordinates[3]))

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξ = T(Iglobal.I[1] - nhalo)
    η = T(Iglobal.I[2] - nhalo)
    ζ = T(Iglobal.I[3] - nhalo)
    x_h[I] = mapping_functions.x1(t, ξ, η, ζ, params)
    y_h[I] = mapping_functions.x2(t, ξ, η, ζ, params)
    z_h[I] = mapping_functions.x3(t, ξ, η, ζ, params)
  end

  copyto!(node_coordinates[1], x_h)
  copyto!(node_coordinates[2], y_h)
  copyto!(node_coordinates[3], z_h)
  return nothing
end

function _compute_unified_centroid_coordinates!(
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  backend,
) where {N}
  _compute_unified_centroid_coordinates!(
    backend, centroid_coordinates, mapping_functions, iterators, nhalo, t, params, Val(N)
  )
end

function _compute_unified_centroid_coordinates!(
  backend::KernelAbstractions.CPU,
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  _compute_unified_centroids_1d_kernel!(backend)(
    centroid_coordinates[1],
    mapping_functions.x1,
    local_domain,
    global_domain,
    nhalo,
    t,
    params;
    ndrange=size(local_domain),
  )
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_centroid_coordinates!(
  backend::KernelAbstractions.CPU,
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  _compute_unified_centroids_2d_kernel!(backend)(
    centroid_coordinates[1],
    centroid_coordinates[2],
    mapping_functions.x1,
    mapping_functions.x2,
    local_domain,
    global_domain,
    nhalo,
    t,
    params;
    ndrange=size(local_domain),
  )
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_centroid_coordinates!(
  backend::KernelAbstractions.CPU,
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  _compute_unified_centroids_3d_kernel!(backend)(
    centroid_coordinates[1],
    centroid_coordinates[2],
    centroid_coordinates[3],
    mapping_functions.x1,
    mapping_functions.x2,
    mapping_functions.x3,
    local_domain,
    global_domain,
    nhalo,
    t,
    params;
    ndrange=size(local_domain),
  )
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_centroid_coordinates!(
  backend::KernelAbstractions.Backend,
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  T = eltype(centroid_coordinates[1])
  half = T(0.5)
  x_h = Array{T}(undef, size(centroid_coordinates[1]))

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξ = T(Iglobal.I[1] - nhalo) + half
    x_h[I] = mapping_functions.x1(t, ξ, params)
  end

  copyto!(centroid_coordinates[1], x_h)
  return nothing
end

function _compute_unified_centroid_coordinates!(
  backend::KernelAbstractions.Backend,
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  T = eltype(centroid_coordinates[1])
  half = T(0.5)
  x_h = Array{T}(undef, size(centroid_coordinates[1]))
  y_h = Array{T}(undef, size(centroid_coordinates[2]))

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξ = T(Iglobal.I[1] - nhalo) + half
    η = T(Iglobal.I[2] - nhalo) + half
    x_h[I] = mapping_functions.x1(t, ξ, η, params)
    y_h[I] = mapping_functions.x2(t, ξ, η, params)
  end

  copyto!(centroid_coordinates[1], x_h)
  copyto!(centroid_coordinates[2], y_h)
  return nothing
end

function _compute_unified_centroid_coordinates!(
  backend::KernelAbstractions.Backend,
  centroid_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  T = eltype(centroid_coordinates[1])
  half = T(0.5)
  x_h = Array{T}(undef, size(centroid_coordinates[1]))
  y_h = Array{T}(undef, size(centroid_coordinates[2]))
  z_h = Array{T}(undef, size(centroid_coordinates[3]))

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξ = T(Iglobal.I[1] - nhalo) + half
    η = T(Iglobal.I[2] - nhalo) + half
    ζ = T(Iglobal.I[3] - nhalo) + half
    x_h[I] = mapping_functions.x1(t, ξ, η, ζ, params)
    y_h[I] = mapping_functions.x2(t, ξ, η, ζ, params)
    z_h[I] = mapping_functions.x3(t, ξ, η, ζ, params)
  end

  copyto!(centroid_coordinates[1], x_h)
  copyto!(centroid_coordinates[2], y_h)
  copyto!(centroid_coordinates[3], z_h)
  return nothing
end

function _compute_unified_face_coordinates!(
  face_coordinates, mapping_functions, iterators, nhalo::Int, t, params, ::Val{N}, backend
) where {N}
  _compute_unified_face_coordinates!(
    backend, face_coordinates, mapping_functions, iterators, nhalo, t, params, Val(N)
  )
end

function _compute_unified_face_coordinates!(
  backend::KernelAbstractions.CPU,
  face_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  _compute_unified_face_coordinates_1d_kernel!(backend)(
    face_coordinates[1][1],
    mapping_functions.x1,
    local_domain,
    global_domain,
    nhalo,
    t,
    params;
    ndrange=size(local_domain),
  )
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_face_coordinates!(
  backend::KernelAbstractions.CPU,
  face_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  for axis in 1:2
    _compute_unified_face_coordinates_2d_kernel!(backend)(
      face_coordinates[axis][1],
      face_coordinates[axis][2],
      mapping_functions.x1,
      mapping_functions.x2,
      local_domain,
      global_domain,
      axis,
      nhalo,
      t,
      params;
      ndrange=size(local_domain),
    )
  end
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_face_coordinates!(
  backend::KernelAbstractions.CPU,
  face_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  for axis in 1:3
    _compute_unified_face_coordinates_3d_kernel!(backend)(
      face_coordinates[axis][1],
      face_coordinates[axis][2],
      face_coordinates[axis][3],
      mapping_functions.x1,
      mapping_functions.x2,
      mapping_functions.x3,
      local_domain,
      global_domain,
      axis,
      nhalo,
      t,
      params;
      ndrange=size(local_domain),
    )
  end
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _compute_unified_face_coordinates!(
  backend::KernelAbstractions.Backend,
  face_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{1},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  T = eltype(face_coordinates[1][1])
  x_h = Array{T}(undef, size(face_coordinates[1][1]))

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξ = T(Iglobal.I[1] - nhalo + 1)
    x_h[I] = mapping_functions.x1(t, ξ, params)
  end

  copyto!(face_coordinates[1][1], x_h)
  return nothing
end

function _compute_unified_face_coordinates!(
  backend::KernelAbstractions.Backend,
  face_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{2},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  T = eltype(face_coordinates[1][1])
  half = T(0.5)

  for axis in 1:2
    x_h = Array{T}(undef, size(face_coordinates[axis][1]))
    y_h = Array{T}(undef, size(face_coordinates[axis][2]))
    @inbounds for I in local_domain
      Iglobal = global_domain[I]
      ξ = T(Iglobal.I[1] - nhalo) + half + (axis == 1 ? half : zero(T))
      η = T(Iglobal.I[2] - nhalo) + half + (axis == 2 ? half : zero(T))
      x_h[I] = mapping_functions.x1(t, ξ, η, params)
      y_h[I] = mapping_functions.x2(t, ξ, η, params)
    end
    copyto!(face_coordinates[axis][1], x_h)
    copyto!(face_coordinates[axis][2], y_h)
  end

  return nothing
end

function _compute_unified_face_coordinates!(
  backend::KernelAbstractions.Backend,
  face_coordinates,
  mapping_functions,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{3},
)
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  T = eltype(face_coordinates[1][1])
  half = T(0.5)

  for axis in 1:3
    x_h = Array{T}(undef, size(face_coordinates[axis][1]))
    y_h = Array{T}(undef, size(face_coordinates[axis][2]))
    z_h = Array{T}(undef, size(face_coordinates[axis][3]))
    @inbounds for I in local_domain
      Iglobal = global_domain[I]
      ξ = T(Iglobal.I[1] - nhalo) + half + (axis == 1 ? half : zero(T))
      η = T(Iglobal.I[2] - nhalo) + half + (axis == 2 ? half : zero(T))
      ζ = T(Iglobal.I[3] - nhalo) + half + (axis == 3 ? half : zero(T))
      x_h[I] = mapping_functions.x1(t, ξ, η, ζ, params)
      y_h[I] = mapping_functions.x2(t, ξ, η, ζ, params)
      z_h[I] = mapping_functions.x3(t, ξ, η, ζ, params)
    end
    copyto!(face_coordinates[axis][1], x_h)
    copyto!(face_coordinates[axis][2], y_h)
    copyto!(face_coordinates[axis][3], z_h)
  end

  return nothing
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{1}, edge, edge_axis::Int, t, ξηζ, params
)
  if edge_axis != 1
    throw(ArgumentError("Invalid 1D face axis: $edge_axis"))
  end

  jinv_edge = edge.Jinv_ᵢ₊½
  norm_edge = edge.norm_Jinv_ᵢ₊½

  jinv_fun = jinv_edge isa NamedTuple ? jinv_edge.ϕᵢ₊½ : jinv_edge
  norm_fun = norm_edge isa NamedTuple ? norm_edge.ϕᵢ₊½ : norm_edge

  G = _as_smatrix(Val(1), jinv_fun(t, ξηζ..., params))
  Ghat = _as_smatrix(Val(1), norm_fun(t, ξηζ..., params))
  Jinv = det(G)
  return G, Ghat, Jinv
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{2}, edge, edge_axis::Int, t, ξηζ, params
)
  if edge_axis == 1
    G = _as_smatrix(Val(2), edge.Jinv_ᵢ₊½(t, ξηζ..., params))
    Ghat = _as_smatrix(Val(2), edge.norm_Jinv_ᵢ₊½(t, ξηζ..., params))
  elseif edge_axis == 2
    G = _as_smatrix(Val(2), edge.Jinv_ⱼ₊½(t, ξηζ..., params))
    Ghat = _as_smatrix(Val(2), edge.norm_Jinv_ⱼ₊½(t, ξηζ..., params))
  else
    throw(ArgumentError("Invalid 2D face axis: $edge_axis"))
  end
  Jinv = det(G)
  return G, Ghat, Jinv
end

@inline function _inverse_and_normalized_edge_metrics(
  ::Val{3}, edge, edge_axis::Int, t, ξηζ, params
)
  if edge_axis == 1
    G = _as_smatrix(Val(3), edge.Jinvᵢ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xᵢ₊½(t, ξηζ..., params) edge.ξ̂yᵢ₊½(t, ξηζ..., params) edge.ξ̂zᵢ₊½(t, ξηζ..., params)
      edge.η̂xᵢ₊½(t, ξηζ..., params) edge.η̂yᵢ₊½(t, ξηζ..., params) edge.η̂zᵢ₊½(t, ξηζ..., params)
      edge.ζ̂xᵢ₊½(t, ξηζ..., params) edge.ζ̂yᵢ₊½(t, ξηζ..., params) edge.ζ̂zᵢ₊½(t, ξηζ..., params)
    ]
  elseif edge_axis == 2
    G = _as_smatrix(Val(3), edge.Jinvⱼ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xⱼ₊½(t, ξηζ..., params) edge.ξ̂yⱼ₊½(t, ξηζ..., params) edge.ξ̂zⱼ₊½(t, ξηζ..., params)
      edge.η̂xⱼ₊½(t, ξηζ..., params) edge.η̂yⱼ₊½(t, ξηζ..., params) edge.η̂zⱼ₊½(t, ξηζ..., params)
      edge.ζ̂xⱼ₊½(t, ξηζ..., params) edge.ζ̂yⱼ₊½(t, ξηζ..., params) edge.ζ̂zⱼ₊½(t, ξηζ..., params)
    ]
  elseif edge_axis == 3
    G = _as_smatrix(Val(3), edge.Jinvₖ₊½(t, ξηζ..., params))
    Ghat = @SMatrix [
      edge.ξ̂xₖ₊½(t, ξηζ..., params) edge.ξ̂yₖ₊½(t, ξηζ..., params) edge.ξ̂zₖ₊½(t, ξηζ..., params)
      edge.η̂xₖ₊½(t, ξηζ..., params) edge.η̂yₖ₊½(t, ξηζ..., params) edge.η̂zₖ₊½(t, ξηζ..., params)
      edge.ζ̂xₖ₊½(t, ξηζ..., params) edge.ζ̂yₖ₊½(t, ξηζ..., params) edge.ζ̂zₖ₊½(t, ξηζ..., params)
    ]
  else
    throw(ArgumentError("Invalid 3D face axis: $edge_axis"))
  end

  Jinv = det(G)
  return G, Ghat, Jinv
end

@kernel function _fill_cell_metric_storage_1d_kernel!(
  forward,
  inverse,
  jacobian,
  inverse_jacobian,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half

  F = _as_smatrix(Val(1), jacobian(t, ξ, params))
  G = _as_smatrix(Val(1), inverse_jacobian(t, ξ, params))
  J = det(F)
  Jinv = det(G)

  forward[I] = Metric(SMatrix{1,1,T,1}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{1,1,T,1}(Tuple(G)), T(J))
end

@kernel function _fill_cell_metric_storage_2d_kernel!(
  forward,
  inverse,
  jacobian,
  inverse_jacobian,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half

  F = _as_smatrix(Val(2), jacobian(t, ξ, η, params))
  G = _as_smatrix(Val(2), inverse_jacobian(t, ξ, η, params))
  J = det(F)
  Jinv = det(G)

  forward[I] = Metric(SMatrix{2,2,T,4}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{2,2,T,4}(Tuple(G)), T(J))
end

@kernel function _fill_cell_metric_storage_3d_kernel!(
  forward,
  inverse,
  jacobian,
  inverse_jacobian,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  F = _as_smatrix(Val(3), jacobian(t, ξ, η, ζ, params))
  G = _as_smatrix(Val(3), inverse_jacobian(t, ξ, η, ζ, params))
  J = det(F)
  Jinv = det(G)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
end

@kernel function _fill_face_metric_storage_1d_axis1_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(Val(1), edge, 1, t, (ξ,), params)
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{1,1,T,1}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{1,1,T,1}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{1,1,T,1}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_2d_axis1_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(Val(2), edge, 1, t, (ξ, η), params)
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{2,2,T,4}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{2,2,T,4}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{2,2,T,4}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_2d_axis2_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(Val(2), edge, 2, t, (ξ, η), params)
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{2,2,T,4}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{2,2,T,4}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{2,2,T,4}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_3d_axis1_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(
    Val(3), edge, 1, t, (ξ, η, ζ), params
  )
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{3,3,T,9}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_3d_axis2_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(
    Val(3), edge, 2, t, (ξ, η, ζ), params
  )
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{3,3,T,9}(Tuple(Ghat)))
end

@kernel function _fill_face_metric_storage_3d_axis3_kernel!(
  forward,
  inverse,
  conserved,
  edge,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  Iglobal = global_domain[I]
  half = T(0.5)
  ξ = Iglobal.I[1] - nhalo + half
  η = Iglobal.I[2] - nhalo + half
  ζ = Iglobal.I[3] - nhalo + half

  G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(
    Val(3), edge, 3, t, (ξ, η, ζ), params
  )
  F = inv(G)
  J = det(F)

  forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
  inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
  conserved[I] = ConservedMetric(SMatrix{3,3,T,9}(Tuple(Ghat)))
end

function _fill_cell_metric_storage!(
  cell_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
  backend,
) where {N,T}
  _fill_cell_metric_storage!(
    backend,
    cell_metric_storage,
    metric_functions_cache,
    iterators,
    nhalo,
    t,
    params,
    Val(N),
    T,
  )
end

function _fill_cell_metric_storage!(
  backend::KernelAbstractions.Backend,
  cell_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full

  forward_h = Array{eltype(cell_metric_storage.forward)}(
    undef, size(cell_metric_storage.forward)
  )
  inverse_h = Array{eltype(cell_metric_storage.inverse)}(
    undef, size(cell_metric_storage.inverse)
  )
  half = T(0.5)

  @inbounds for I in local_domain
    Iglobal = global_domain[I]
    ξηζ = ntuple(d -> T(Iglobal.I[d] - nhalo) + half, N)

    F = _as_smatrix(Val(N), metric_functions_cache.forward.jacobian(t, ξηζ..., params))
    G = _as_smatrix(Val(N), metric_functions_cache.inverse.Jinv(t, ξηζ..., params))
    J = det(F)
    Jinv = det(G)

    forward_h[I] = Metric(SMatrix{N,N,T,N * N}(Tuple(F)), T(J))
    inverse_h[I] = Metric(SMatrix{N,N,T,N * N}(Tuple(G)), T(J))
  end

  copyto!(cell_metric_storage.forward, forward_h)
  copyto!(cell_metric_storage.inverse, inverse_h)
  return nothing
end

function _fill_cell_metric_storage!(
  backend::KernelAbstractions.CPU,
  cell_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full

  if N == 1
    _fill_cell_metric_storage_1d_kernel!(backend)(
      cell_metric_storage.forward,
      cell_metric_storage.inverse,
      metric_functions_cache.forward.jacobian,
      metric_functions_cache.inverse.Jinv,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 2
    _fill_cell_metric_storage_2d_kernel!(backend)(
      cell_metric_storage.forward,
      cell_metric_storage.inverse,
      metric_functions_cache.forward.jacobian,
      metric_functions_cache.inverse.Jinv,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 3
    _fill_cell_metric_storage_3d_kernel!(backend)(
      cell_metric_storage.forward,
      cell_metric_storage.inverse,
      metric_functions_cache.forward.jacobian,
      metric_functions_cache.inverse.Jinv,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  else
    throw(ArgumentError("Unsupported metric dimension N=$N"))
  end
  KernelAbstractions.synchronize(backend)
  return nothing
end

function _fill_face_metric_storage!(
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
  backend,
) where {N,T}
  _fill_face_metric_storage!(
    backend,
    face_metric_storage,
    metric_functions_cache,
    iterators,
    nhalo,
    t,
    params,
    Val(N),
    T,
  )
  return nothing
end

function _fill_face_metric_storage!(
  backend::KernelAbstractions.Backend,
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  edge = metric_functions_cache.edge
  half = T(0.5)

  for axis in 1:N
    forward_h = Array{eltype(face_metric_storage[axis].forward)}(
      undef, size(face_metric_storage[axis].forward)
    )
    inverse_h = Array{eltype(face_metric_storage[axis].inverse)}(
      undef, size(face_metric_storage[axis].inverse)
    )
    conserved_h = Array{eltype(face_metric_storage[axis].conserved)}(
      undef, size(face_metric_storage[axis].conserved)
    )

    @inbounds for I in local_domain
      Iglobal = global_domain[I]
      ξηζ = ntuple(d -> T(Iglobal.I[d] - nhalo) + half, N)

      G, Ghat, Jinv = _inverse_and_normalized_edge_metrics(
        Val(N), edge, axis, t, ξηζ, params
      )
      F = inv(G)
      J = det(F)

      forward_h[I] = Metric(SMatrix{N,N,T,N * N}(Tuple(F)), T(J))
      inverse_h[I] = Metric(SMatrix{N,N,T,N * N}(Tuple(G)), T(J))
      conserved_h[I] = ConservedMetric(SMatrix{N,N,T,N * N}(Tuple(Ghat)))
    end

    copyto!(face_metric_storage[axis].forward, forward_h)
    copyto!(face_metric_storage[axis].inverse, inverse_h)
    copyto!(face_metric_storage[axis].conserved, conserved_h)
  end
  return nothing
end

@inline function _padded_global_cell_index(
  global_domain, i::Int, j::Int, k::Int, ni::Int, nj::Int, nk::Int
)
  ic = clamp(i, 1, ni)
  jc = clamp(j, 1, nj)
  kc = clamp(k, 1, nk)
  Iglobal = global_domain[CartesianIndex(ic, jc, kc)]
  return (
    Iglobal.I[1] + i - ic,
    Iglobal.I[2] + j - jc,
    Iglobal.I[3] + k - kc,
  )
end

@inline function _shift_array_axis(
  i::Int, j::Int, k::Int, axis::Int, offset::Int
)
  if axis == 1
    return i + offset, j, k
  elseif axis == 2
    return i, j + offset, k
  elseif axis == 3
    return i, j, k + offset
  else
    throw(ArgumentError("Invalid 3D derivative axis: $axis"))
  end
end

@inline _shift_array_axis(i::Int, j::Int, k::Int, ::Val{1}, offset::Int) =
  (i + offset, j, k)
@inline _shift_array_axis(i::Int, j::Int, k::Int, ::Val{2}, offset::Int) =
  (i, j + offset, k)
@inline _shift_array_axis(i::Int, j::Int, k::Int, ::Val{3}, offset::Int) =
  (i, j, k + offset)

@inline function _ad_thomas_lombard_potential(edge, potential_id::Int, t, ξ, η, ζ, params)
  F = edge.jacobian(t, ξ, η, ζ, params)

  return _ad_thomas_lombard_potential_from_values(
    F,
    edge.x(t, ξ, η, ζ, params),
    edge.y(t, ξ, η, ζ, params),
    edge.z(t, ξ, η, ζ, params),
    potential_id,
  )
end

@inline function _ad_thomas_lombard_potential_from_values(F, x, y, z, potential_id::Int)
  if potential_id == 1
    return F[2, 2] * z
  elseif potential_id == 2
    return F[2, 3] * z
  elseif potential_id == 3
    return F[2, 1] * z
  elseif potential_id == 4
    return F[3, 2] * x
  elseif potential_id == 5
    return F[3, 3] * x
  elseif potential_id == 6
    return F[3, 1] * x
  elseif potential_id == 7
    return F[1, 2] * y
  elseif potential_id == 8
    return F[1, 3] * y
  elseif potential_id == 9
    return F[1, 1] * y
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard potential id: $potential_id"))
  end
end

@inline function _ad_thomas_lombard_potential_vector(
  edge, potential_ids::NTuple{3,Int}, t, ξ, η, ζ, params
)
  F = edge.jacobian(t, ξ, η, ζ, params)
  x = edge.x(t, ξ, η, ζ, params)
  y = edge.y(t, ξ, η, ζ, params)
  z = edge.z(t, ξ, η, ζ, params)
  return SVector{3}(
    ntuple(
      n -> _ad_thomas_lombard_potential_from_values(F, x, y, z, potential_ids[n]),
      Val(3),
    ),
  )
end

@inline function _ad_thomas_lombard_potential_vector(
  edge, potential_ids::Val, t, ξ, η, ζ, params
)
  return _ad_thomas_lombard_potential_vector_from_axis_derivative(
    edge, potential_ids, t, ξ, η, ζ, params
  )
end

@generated function _ad_thomas_lombard_seed_axis(::Val{Axis}, ξ, η, ζ) where {Axis}
  Axis in (1, 2, 3) || error("Invalid AD Thomas-Lombard seed axis: $Axis")
  vars = (:ξ, :η, :ζ)
  coord = vars[Axis]
  seeded = collect(vars)
  seeded[Axis] = :seeded_coord
  return quote
    coord_value = $coord + zero(ξ) + zero(η) + zero(ζ)
    tag = ForwardDiff.Tag(_ad_thomas_lombard_seed_axis, typeof(coord_value))
    tag_type = typeof(tag)
    seeded_coord = ForwardDiff.Dual{tag_type}(coord_value, oneunit(coord_value))
    return ($(seeded[1]), $(seeded[2]), $(seeded[3]), tag_type)
  end
end

@inline _ad_thomas_lombard_pair_axis_a_tag() = nothing
@inline _ad_thomas_lombard_pair_axis_b_tag() = nothing

@inline function _ad_thomas_lombard_seed_second_axis(
  ::Val{1}, tag_function, ξ, η, ζ
)
  coord_value = ξ
  tag1 = ForwardDiff.Tag(tag_function, typeof(coord_value))
  tag1_type = typeof(tag1)
  coord1 = ForwardDiff.Dual{tag1_type}(coord_value, oneunit(coord_value))
  tag2 = ForwardDiff.Tag(tag_function, typeof(coord1))
  tag2_type = typeof(tag2)
  coord2 = ForwardDiff.Dual{tag2_type}(coord1, oneunit(coord1))
  return coord2, η, ζ, tag1_type, tag2_type
end

@inline function _ad_thomas_lombard_seed_second_axis(
  ::Val{2}, tag_function, ξ, η, ζ
)
  coord_value = η
  tag1 = ForwardDiff.Tag(tag_function, typeof(coord_value))
  tag1_type = typeof(tag1)
  coord1 = ForwardDiff.Dual{tag1_type}(coord_value, oneunit(coord_value))
  tag2 = ForwardDiff.Tag(tag_function, typeof(coord1))
  tag2_type = typeof(tag2)
  coord2 = ForwardDiff.Dual{tag2_type}(coord1, oneunit(coord1))
  return ξ, coord2, ζ, tag1_type, tag2_type
end

@inline function _ad_thomas_lombard_seed_second_axis(
  ::Val{3}, tag_function, ξ, η, ζ
)
  coord_value = ζ
  tag1 = ForwardDiff.Tag(tag_function, typeof(coord_value))
  tag1_type = typeof(tag1)
  coord1 = ForwardDiff.Dual{tag1_type}(coord_value, oneunit(coord_value))
  tag2 = ForwardDiff.Tag(tag_function, typeof(coord1))
  tag2_type = typeof(tag2)
  coord2 = ForwardDiff.Dual{tag2_type}(coord1, oneunit(coord1))
  return ξ, η, coord2, tag1_type, tag2_type
end

@inline function _ad_thomas_lombard_coordinate_axis_values(
  edge, axis::Val, t, ξ, η, ζ, params
)
  ξd, ηd, ζd, tag_type = _ad_thomas_lombard_seed_axis(axis, ξ, η, ζ)
  xd = edge.x(t, ξd, ηd, ζd, params)
  yd = edge.y(t, ξd, ηd, ζd, params)
  zd = edge.z(t, ξd, ηd, ζd, params)
  x = _dual_value_for_tag(tag_type, xd)
  y = _dual_value_for_tag(tag_type, yd)
  z = _dual_value_for_tag(tag_type, zd)
  dx = _dual_derivative_for_tag(tag_type, xd)
  dy = _dual_derivative_for_tag(tag_type, yd)
  dz = _dual_derivative_for_tag(tag_type, zd)
  return x, y, z, dx, dy, dz
end

@inline function _ad_thomas_lombard_potential_vector_from_axis_derivative(
  edge, ::Val{(1, 4, 7)}, t, ξ, η, ζ, params
)
  x, y, z, dx, dy, dz = _ad_thomas_lombard_coordinate_axis_values(
    edge, Val(2), t, ξ, η, ζ, params
  )
  return SVector(dy * z, dz * x, dx * y)
end

@inline function _ad_thomas_lombard_potential_vector_from_axis_derivative(
  edge, ::Val{(2, 5, 8)}, t, ξ, η, ζ, params
)
  x, y, z, dx, dy, dz = _ad_thomas_lombard_coordinate_axis_values(
    edge, Val(3), t, ξ, η, ζ, params
  )
  return SVector(dy * z, dz * x, dx * y)
end

@inline function _ad_thomas_lombard_potential_vector_from_axis_derivative(
  edge, ::Val{(3, 6, 9)}, t, ξ, η, ζ, params
)
  x, y, z, dx, dy, dz = _ad_thomas_lombard_coordinate_axis_values(
    edge, Val(1), t, ξ, η, ζ, params
  )
  return SVector(dy * z, dz * x, dx * y)
end

@inline _dual_value_for_tag(::Type{Tag}, y) where {Tag} = ForwardDiff.value(Tag, y)
@inline _dual_derivative_for_tag(::Type{Tag}, y) where {Tag} =
  ForwardDiff.partials(Tag, y, 1)

@inline function _dual_value_for_tag(::Type{Tag}, y::SVector{N}) where {Tag,N}
  return SVector{N}(ntuple(i -> _dual_value_for_tag(Tag, y[i]), Val(N)))
end

@inline function _dual_derivative_for_tag(::Type{Tag}, y::SVector{N}) where {Tag,N}
  return SVector{N}(ntuple(i -> _dual_derivative_for_tag(Tag, y[i]), Val(N)))
end

@inline function _dual_value_for_tag(::Type{Tag}, y::Tuple) where {Tag}
  return map(v -> _dual_value_for_tag(Tag, v), y)
end

@inline function _dual_derivative_for_tag(::Type{Tag}, y::Tuple) where {Tag}
  return map(v -> _dual_derivative_for_tag(Tag, v), y)
end

@inline function _forwarddiff_extract_value_derivative_and_second_derivative(
  ::Type{T1}, ::Type{T2}, ydual2
) where {T1,T2}
  ydual1 = _dual_value_for_tag(T2, ydual2)
  value = _dual_value_for_tag(T1, ydual1)
  deriv = _dual_derivative_for_tag(T1, ydual1)
  second_deriv = _dual_derivative_for_tag(T1, _dual_derivative_for_tag(T2, ydual2))
  return value, deriv, second_deriv
end

@inline function _forwarddiff_value_derivative_and_second_derivative(f, x)
  tag1 = ForwardDiff.Tag(f, typeof(x))
  T1 = typeof(tag1)
  xdual1 = ForwardDiff.Dual{T1}(x, oneunit(x))
  tag2 = ForwardDiff.Tag(f, typeof(xdual1))
  T2 = typeof(tag2)
  xdual2 = ForwardDiff.Dual{T2}(xdual1, oneunit(xdual1))
  ydual2 = f(xdual2)
  return _forwarddiff_extract_value_derivative_and_second_derivative(T1, T2, ydual2)
end

function _ad_thomas_lombard_potential_axis_derivs(
  edge, potential_id::Int, axis::Int, t, ξ, η, ζ, params
)
  backend = edge.diff_backend
  if axis == 1
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential(edge, potential_id, t, s, η, ζ, params),
      backend,
      ξ,
    )
  elseif axis == 2
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential(edge, potential_id, t, ξ, s, ζ, params),
      backend,
      η,
    )
  elseif axis == 3
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential(edge, potential_id, t, ξ, η, s, params),
      backend,
      ζ,
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard potential derivative axis: $axis"))
  end
end

function _ad_thomas_lombard_potential_axis_derivs_vector(
  edge, potential_ids, axis::Int, t, ξ, η, ζ, params
)
  return _ad_thomas_lombard_potential_axis_derivs_vector(
    edge, potential_ids, Val(axis), t, ξ, η, ζ, params
  )
end

function _ad_thomas_lombard_potential_axis_derivs_vector(
  edge, potential_ids, ::Val{Axis}, t, ξ, η, ζ, params
) where {Axis}
  if Axis == 1
    return _forwarddiff_value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential_vector(edge, potential_ids, t, s, η, ζ, params),
      ξ,
    )
  elseif Axis == 2
    return _forwarddiff_value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential_vector(edge, potential_ids, t, ξ, s, ζ, params),
      η,
    )
  elseif Axis == 3
    return _forwarddiff_value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_potential_vector(edge, potential_ids, t, ξ, η, s, params),
      ζ,
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard potential derivative axis: $Axis"))
  end
end

function _ad_thomas_lombard_edge_potential(
  edge, potential_id::Int, edge_axis::Int, t, ξ, η, ζ, params
)
  ϕᵢ, ϕaᵢ, ϕaaᵢ = _ad_thomas_lombard_potential_axis_derivs(
    edge, potential_id, edge_axis, t, ξ, η, ζ, params
  )

  if edge_axis == 1
    ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁ = _ad_thomas_lombard_potential_axis_derivs(
      edge, potential_id, edge_axis, t, ξ + one(ξ), η, ζ, params
    )
  elseif edge_axis == 2
    ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁ = _ad_thomas_lombard_potential_axis_derivs(
      edge, potential_id, edge_axis, t, ξ, η + one(η), ζ, params
    )
  elseif edge_axis == 3
    ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁ = _ad_thomas_lombard_potential_axis_derivs(
      edge, potential_id, edge_axis, t, ξ, η, ζ + one(ζ), params
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard edge axis: $edge_axis"))
  end

  return _edge_reconstruct(
    ϕᵢ, ϕaᵢ, ϕaaᵢ, ϕᵢ₊₁, ϕaᵢ₊₁, ϕaaᵢ₊₁, EdgeInterpolationOrder3()
  )
end

function _ad_thomas_lombard_edge_potential_derivs(
  edge, potential_id::Int, edge_axis::Int, deriv_axis::Int, t, ξ, η, ζ, params
)
  backend = edge.diff_backend
  if deriv_axis == 1
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_edge_potential(
        edge, potential_id, edge_axis, t, s, η, ζ, params
      ),
      backend,
      ξ,
    )
  elseif deriv_axis == 2
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_edge_potential(
        edge, potential_id, edge_axis, t, ξ, s, ζ, params
      ),
      backend,
      η,
    )
  elseif deriv_axis == 3
    return value_derivative_and_second_derivative(
      s -> _ad_thomas_lombard_edge_potential(
        edge, potential_id, edge_axis, t, ξ, η, s, params
      ),
      backend,
      ζ,
    )
  else
    throw(ArgumentError("Invalid AD Thomas-Lombard derivative axis: $deriv_axis"))
  end
end

function _ad_thomas_lombard_potential_pair_derivs_vector(
  edge, potential_ids, axis_a::Int, axis_b::Int, t, ξ, η, ζ, params
)
  return _ad_thomas_lombard_potential_pair_derivs_vector(
    edge, potential_ids, Val(axis_a), Val(axis_b), t, ξ, η, ζ, params
  )
end

function _ad_thomas_lombard_potential_pair_derivs_vector(
  edge, potential_ids, ::Val{AxisA}, ::Val{AxisB}, t, ξ, η, ζ, params
) where {AxisA,AxisB}
  ξb, ηb, ζb, Tb1, Tb2 = _ad_thomas_lombard_seed_second_axis(
    Val(AxisB), _ad_thomas_lombard_pair_axis_b_tag, ξ, η, ζ
  )
  ξab, ηab, ζab, Ta1, Ta2 = _ad_thomas_lombard_seed_second_axis(
    Val(AxisA), _ad_thomas_lombard_pair_axis_a_tag, ξb, ηb, ζb
  )

  Pdual = _ad_thomas_lombard_potential_vector(
    edge, potential_ids, t, ξab, ηab, ζab, params
  )
  P_bdual, Pa_bdual, Paa_bdual =
    _forwarddiff_extract_value_derivative_and_second_derivative(Ta1, Ta2, Pdual)

  P, Pb, Pbb =
    _forwarddiff_extract_value_derivative_and_second_derivative(Tb1, Tb2, P_bdual)
  Pa, Pab, Pabb =
    _forwarddiff_extract_value_derivative_and_second_derivative(Tb1, Tb2, Pa_bdual)
  Paa, Paab, Paabb =
    _forwarddiff_extract_value_derivative_and_second_derivative(Tb1, Tb2, Paa_bdual)
  return P, Pa, Paa, Pb, Pab, Paab, Pbb, Pabb, Paabb
end

function _fill_ad_thomas_lombard_edge_potential_derivs_3d!(
  Q,
  Qb,
  Qbb,
  edge,
  potential_id::Int,
  edge_axis::Int,
  deriv_axis::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  ni, nj, nk = size(local_domain)
  half = T(0.5)
  irange = deriv_axis == 1 ? (1:(ni + 2)) : (2:(ni + 1))
  jrange = deriv_axis == 2 ? (1:(nj + 2)) : (2:(nj + 1))
  krange = deriv_axis == 3 ? (1:(nk + 2)) : (2:(nk + 1))

  Threads.@threads for kp in krange
    @inbounds for jp in jrange, ip in irange
      i = ip - 1
      j = jp - 1
      k = kp - 1
      Iglobal = _padded_global_cell_index(global_domain, i, j, k, ni, nj, nk)
      ξ = T(Iglobal[1] - nhalo) + half
      η = T(Iglobal[2] - nhalo) + half
      ζ = T(Iglobal[3] - nhalo) + half

      q, qb, qbb = _ad_thomas_lombard_edge_potential_derivs(
        edge, potential_id, edge_axis, deriv_axis, t, ξ, η, ζ, params
      )
      Q[ip, jp, kp] = T(q)
      Qb[ip, jp, kp] = T(qb)
      Qbb[ip, jp, kp] = T(qbb)
    end
  end

  return nothing
end

function _fill_ad_thomas_lombard_potential_pair_table_3d!(
  table,
  edge,
  potential_ids,
  axis_a::Int,
  axis_b::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  return _fill_ad_thomas_lombard_potential_pair_table_3d!(
    table,
    edge,
    potential_ids,
    Val(axis_a),
    Val(axis_b),
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )
end

function _fill_ad_thomas_lombard_potential_pair_table_3d!(
  table,
  edge,
  potential_ids,
  ::Val{AxisA},
  ::Val{AxisB},
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T,AxisA,AxisB}
  ni, nj, nk = size(local_domain)
  half = T(0.5)
  irange = AxisA == 1 || AxisB == 1 ? (1:(ni + 2)) : (2:(ni + 1))
  jrange = AxisA == 2 || AxisB == 2 ? (1:(nj + 2)) : (2:(nj + 1))
  krange = AxisA == 3 || AxisB == 3 ? (1:(nk + 2)) : (2:(nk + 1))

  Threads.@threads for kp in krange
    @inbounds for jp in jrange, ip in irange
      i = ip - 1
      j = jp - 1
      k = kp - 1
      Iglobal = _padded_global_cell_index(global_domain, i, j, k, ni, nj, nk)
      ξ = T(Iglobal[1] - nhalo) + half
      η = T(Iglobal[2] - nhalo) + half
      ζ = T(Iglobal[3] - nhalo) + half

      vals = _ad_thomas_lombard_potential_pair_derivs_vector(
        edge, potential_ids, Val(AxisA), Val(AxisB), t, ξ, η, ζ, params
      )
      table[1][ip, jp, kp] = SVector{3,T}(vals[1])
      table[2][ip, jp, kp] = SVector{3,T}(vals[2])
      table[3][ip, jp, kp] = SVector{3,T}(vals[3])
      table[4][ip, jp, kp] = SVector{3,T}(vals[4])
      table[5][ip, jp, kp] = SVector{3,T}(vals[5])
      table[6][ip, jp, kp] = SVector{3,T}(vals[6])
      table[7][ip, jp, kp] = SVector{3,T}(vals[7])
      table[8][ip, jp, kp] = SVector{3,T}(vals[8])
      table[9][ip, jp, kp] = SVector{3,T}(vals[9])
    end
  end

  return nothing
end

@inline function _array_derivative_from_ad_thomas_lombard_edge_potential(
  Q, Qb, Qbb, deriv_axis::Int, i::Int, j::Int, k::Int
)
  ip, jp, kp = _shift_array_axis(i, j, k, deriv_axis, 1)
  im, jm, km = _shift_array_axis(i, j, k, deriv_axis, -1)

  edge_high = _edge_reconstruct(
    Q[i, j, k],
    Qb[i, j, k],
    Qbb[i, j, k],
    Q[ip, jp, kp],
    Qb[ip, jp, kp],
    Qbb[ip, jp, kp],
    EdgeInterpolationOrder3(),
  )
  edge_low = _edge_reconstruct(
    Q[im, jm, km],
    Qb[im, jm, km],
    Qbb[im, jm, km],
    Q[i, j, k],
    Qb[i, j, k],
    Qbb[i, j, k],
    EdgeInterpolationOrder3(),
  )

  return edge_high - edge_low
end

@inline function _edge_reconstruct_from_pair_table(
  table,
  value_id::Int,
  normal_deriv_id::Int,
  normal_second_deriv_id::Int,
  normal_axis::Int,
  i::Int,
  j::Int,
  k::Int,
)
  return _edge_reconstruct_from_pair_table(
    table,
    value_id,
    normal_deriv_id,
    normal_second_deriv_id,
    Val(normal_axis),
    i,
    j,
    k,
  )
end

@inline function _edge_reconstruct_from_pair_table(
  table,
  value_id::Int,
  normal_deriv_id::Int,
  normal_second_deriv_id::Int,
  normal_axis::Val,
  i::Int,
  j::Int,
  k::Int,
)
  ip, jp, kp = _shift_array_axis(i, j, k, normal_axis, 1)
  return _edge_reconstruct(
    table[value_id][i, j, k],
    table[normal_deriv_id][i, j, k],
    table[normal_second_deriv_id][i, j, k],
    table[value_id][ip, jp, kp],
    table[normal_deriv_id][ip, jp, kp],
    table[normal_second_deriv_id][ip, jp, kp],
    EdgeInterpolationOrder3(),
  )
end

@inline function _edge_potential_values_from_pair_table(
  table,
  normal_axis::Int,
  normal_is_axis_a::Bool,
  i::Int,
  j::Int,
  k::Int,
)
  return _edge_potential_values_from_pair_table(
    table, Val(normal_axis), Val(normal_is_axis_a), i, j, k
  )
end

@inline function _edge_potential_values_from_pair_table(
  table,
  normal_axis::Val,
  ::Val{NormalIsAxisA},
  i::Int,
  j::Int,
  k::Int,
) where {NormalIsAxisA}
  if NormalIsAxisA
    Q = _edge_reconstruct_from_pair_table(table, 1, 2, 3, normal_axis, i, j, k)
    Qd = _edge_reconstruct_from_pair_table(table, 4, 5, 6, normal_axis, i, j, k)
    Qdd = _edge_reconstruct_from_pair_table(table, 7, 8, 9, normal_axis, i, j, k)
  else
    Q = _edge_reconstruct_from_pair_table(table, 1, 4, 7, normal_axis, i, j, k)
    Qd = _edge_reconstruct_from_pair_table(table, 2, 5, 8, normal_axis, i, j, k)
    Qdd = _edge_reconstruct_from_pair_table(table, 3, 6, 9, normal_axis, i, j, k)
  end
  return Q, Qd, Qdd
end

@inline function _array_derivative_from_pair_table(
  table,
  normal_axis::Int,
  deriv_axis::Int,
  normal_is_axis_a::Bool,
  i::Int,
  j::Int,
  k::Int,
)
  return _array_derivative_from_pair_table(
    table,
    Val(normal_axis),
    Val(deriv_axis),
    Val(normal_is_axis_a),
    i,
    j,
    k,
  )
end

@inline function _array_derivative_from_pair_table(
  table,
  normal_axis::Val,
  deriv_axis::Val,
  normal_is_axis_a::Val,
  i::Int,
  j::Int,
  k::Int,
)
  ip, jp, kp = _shift_array_axis(i, j, k, deriv_axis, 1)
  im, jm, km = _shift_array_axis(i, j, k, deriv_axis, -1)

  Q, Qd, Qdd = _edge_potential_values_from_pair_table(
    table, normal_axis, normal_is_axis_a, i, j, k
  )
  Qp, Qdp, Qddp = _edge_potential_values_from_pair_table(
    table, normal_axis, normal_is_axis_a, ip, jp, kp
  )
  Qm, Qdm, Qddm = _edge_potential_values_from_pair_table(
    table, normal_axis, normal_is_axis_a, im, jm, km
  )

  edge_high = _edge_reconstruct(Q, Qd, Qdd, Qp, Qdp, Qddp, EdgeInterpolationOrder3())
  edge_low = _edge_reconstruct(Qm, Qdm, Qddm, Q, Qd, Qdd, EdgeInterpolationOrder3())
  return edge_high - edge_low
end

function _accumulate_ad_thomas_lombard_term_3d!(
  out,
  Q,
  Qb,
  Qbb,
  sign,
  edge,
  potential_id::Int,
  edge_axis::Int,
  deriv_axis::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  _fill_ad_thomas_lombard_edge_potential_derivs_3d!(
    Q,
    Qb,
    Qbb,
    edge,
    potential_id,
    edge_axis,
    deriv_axis,
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )

  Threads.@threads for I in local_domain
    @inbounds begin
      i, j, k = I.I
      out[I] += sign * _array_derivative_from_ad_thomas_lombard_edge_potential(
        Q, Qb, Qbb, deriv_axis, i + 1, j + 1, k + 1
      )
    end
  end

  return nothing
end

function _accumulate_ad_thomas_lombard_pair_orientations_3d!(
  out_a,
  sign_a,
  out_b,
  sign_b,
  table,
  axis_a::Val,
  axis_b::Val,
  local_domain,
)
  Threads.@threads for I in local_domain
    @inbounds begin
      i, j, k = I.I
      ip = i + 1
      jp = j + 1
      kp = k + 1
      dQa = _array_derivative_from_pair_table(
        table, axis_a, axis_b, Val(true), ip, jp, kp
      )
      dQb = _array_derivative_from_pair_table(
        table, axis_b, axis_a, Val(false), ip, jp, kp
      )
      out_a[I] += sign_a * dQa
      out_b[I] += sign_b * dQb
    end
  end

  return nothing
end

function _accumulate_ad_thomas_lombard_pair_3d!(
  out_a,
  sign_a,
  out_b,
  sign_b,
  table,
  edge,
  potential_ids,
  axis_a::Int,
  axis_b::Int,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  return _accumulate_ad_thomas_lombard_pair_3d!(
    out_a,
    sign_a,
    out_b,
    sign_b,
    table,
    edge,
    potential_ids,
    Val(axis_a),
    Val(axis_b),
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )
end

function _accumulate_ad_thomas_lombard_pair_3d!(
  out_a,
  sign_a,
  out_b,
  sign_b,
  table,
  edge,
  potential_ids,
  axis_a::Val,
  axis_b::Val,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  _fill_ad_thomas_lombard_potential_pair_table_3d!(
    table,
    edge,
    potential_ids,
    axis_a,
    axis_b,
    local_domain,
    global_domain,
    nhalo,
    t,
    params,
    T,
  )
  _accumulate_ad_thomas_lombard_pair_orientations_3d!(
    out_a,
    sign_a,
    out_b,
    sign_b,
    table,
    axis_a,
    axis_b,
    local_domain,
  )

  return nothing
end

function _fill_face_metric_storage_ad_thomas_lombard_3d_cpu!(
  face_metric_storage,
  metric_functions_cache,
  local_domain,
  global_domain,
  nhalo::Int,
  t,
  params,
  ::Type{T},
) where {T}
  edge = metric_functions_cache.edge
  zeroT = zero(T)
  half = T(0.5)
  row_size = size(local_domain)

  ni, nj, nk = size(local_domain)
  pad_size = (ni + 2, nj + 2, nk + 2)
  if edge.diff_backend isa AutoForwardDiff
    ξrow = fill(zero(SVector{3,T}), row_size)
    ηrow = fill(zero(SVector{3,T}), row_size)
    ζrow = fill(zero(SVector{3,T}), row_size)
    table = ntuple(_ -> Array{SVector{3,T}}(undef, pad_size), Val(9))
    _accumulate_ad_thomas_lombard_pair_3d!(
      ξrow,
      one(T),
      ζrow,
      -one(T),
      table,
      edge,
      Val((1, 4, 7)),
      Val(1),
      Val(3),
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )
    _accumulate_ad_thomas_lombard_pair_3d!(
      ξrow,
      -one(T),
      ηrow,
      one(T),
      table,
      edge,
      Val((2, 5, 8)),
      Val(1),
      Val(2),
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )
    _accumulate_ad_thomas_lombard_pair_3d!(
      ηrow,
      -one(T),
      ζrow,
      one(T),
      table,
      edge,
      Val((3, 6, 9)),
      Val(2),
      Val(3),
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )

    Threads.@threads for I in local_domain
      @inbounds begin
        Iglobal = global_domain[I]
        ξc = T(Iglobal.I[1] - nhalo) + half
        ηc = T(Iglobal.I[2] - nhalo) + half
        ζc = T(Iglobal.I[3] - nhalo) + half

        for axis in 1:3
          ξf, ηf, ζf = _face_center_3d(axis, ξc, ηc, ζc)
          F = _as_smatrix(Val(3), edge.jacobian(t, ξf, ηf, ζf, params))
          G = inv(F)
          J = det(F)

          face_metric_storage[axis].forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
          face_metric_storage[axis].inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
        end

        ξ = ξrow[I]
        η = ηrow[I]
        ζ = ζrow[I]
        face_metric_storage[1].conserved[I] = ConservedMetric(
          @SMatrix [ξ[1] ξ[2] ξ[3]; zeroT zeroT zeroT; zeroT zeroT zeroT]
        )
        face_metric_storage[2].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; η[1] η[2] η[3]; zeroT zeroT zeroT]
        )
        face_metric_storage[3].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; zeroT zeroT zeroT; ζ[1] ζ[2] ζ[3]]
        )
      end
    end
  else
    ξx = fill(zeroT, row_size)
    ξy = fill(zeroT, row_size)
    ξz = fill(zeroT, row_size)
    ηx = fill(zeroT, row_size)
    ηy = fill(zeroT, row_size)
    ηz = fill(zeroT, row_size)
    ζx = fill(zeroT, row_size)
    ζy = fill(zeroT, row_size)
    ζz = fill(zeroT, row_size)
    Q = Array{T}(undef, pad_size)
    Qb = Array{T}(undef, pad_size)
    Qbb = Array{T}(undef, pad_size)
    terms = (
      (ξx, one(T), 1, 1, 3),
      (ξx, -one(T), 2, 1, 2),
      (ηx, one(T), 2, 2, 1),
      (ηx, -one(T), 3, 2, 3),
      (ζx, one(T), 3, 3, 2),
      (ζx, -one(T), 1, 3, 1),
      (ξy, one(T), 4, 1, 3),
      (ξy, -one(T), 5, 1, 2),
      (ηy, one(T), 5, 2, 1),
      (ηy, -one(T), 6, 2, 3),
      (ζy, one(T), 6, 3, 2),
      (ζy, -one(T), 4, 3, 1),
      (ξz, one(T), 7, 1, 3),
      (ξz, -one(T), 8, 1, 2),
      (ηz, one(T), 8, 2, 1),
      (ηz, -one(T), 9, 2, 3),
      (ζz, one(T), 9, 3, 2),
      (ζz, -one(T), 7, 3, 1),
    )
    for (out, sign, potential_id, edge_axis, deriv_axis) in terms
      _accumulate_ad_thomas_lombard_term_3d!(
        out,
        Q,
        Qb,
        Qbb,
        sign,
        edge,
        potential_id,
        edge_axis,
        deriv_axis,
        local_domain,
        global_domain,
        nhalo,
        t,
        params,
        T,
      )
    end

    Threads.@threads for I in local_domain
      @inbounds begin
        Iglobal = global_domain[I]
        ξc = T(Iglobal.I[1] - nhalo) + half
        ηc = T(Iglobal.I[2] - nhalo) + half
        ζc = T(Iglobal.I[3] - nhalo) + half

        for axis in 1:3
          ξf, ηf, ζf = _face_center_3d(axis, ξc, ηc, ζc)
          F = _as_smatrix(Val(3), edge.jacobian(t, ξf, ηf, ζf, params))
          G = inv(F)
          J = det(F)

          face_metric_storage[axis].forward[I] = Metric(SMatrix{3,3,T,9}(Tuple(F)), T(J))
          face_metric_storage[axis].inverse[I] = Metric(SMatrix{3,3,T,9}(Tuple(G)), T(J))
        end

        face_metric_storage[1].conserved[I] = ConservedMetric(
          @SMatrix [ξx[I] ξy[I] ξz[I]; zeroT zeroT zeroT; zeroT zeroT zeroT]
        )
        face_metric_storage[2].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; ηx[I] ηy[I] ηz[I]; zeroT zeroT zeroT]
        )
        face_metric_storage[3].conserved[I] = ConservedMetric(
          @SMatrix [zeroT zeroT zeroT; zeroT zeroT zeroT; ζx[I] ζy[I] ζz[I]]
        )
      end
    end
  end

  return nothing
end

@inline function _face_center_3d(face_axis::Int, ξc::T, ηc::T, ζc::T) where {T}
  half = T(0.5)
  if face_axis == 1
    return ξc + half, ηc, ζc
  elseif face_axis == 2
    return ξc, ηc + half, ζc
  elseif face_axis == 3
    return ξc, ηc, ζc + half
  else
    throw(ArgumentError("Invalid 3D face axis: $face_axis"))
  end
end


function _fill_face_metric_storage!(
  backend::KernelAbstractions.CPU,
  face_metric_storage,
  metric_functions_cache,
  iterators,
  nhalo::Int,
  t,
  params,
  ::Val{N},
  ::Type{T},
) where {N,T}
  local_domain = iterators.cell.full
  global_domain = iterators.global_domain.cell.full
  edge = metric_functions_cache.edge

  if N == 3 && hasproperty(edge, :ad_thomas_lombard_metric)
    _fill_face_metric_storage_ad_thomas_lombard_3d_cpu!(
      face_metric_storage,
      metric_functions_cache,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T,
    )
    return nothing
  end

  if N == 1
    _fill_face_metric_storage_1d_axis1_kernel!(backend)(
      face_metric_storage[1].forward,
      face_metric_storage[1].inverse,
      face_metric_storage[1].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 2
    _fill_face_metric_storage_2d_axis1_kernel!(backend)(
      face_metric_storage[1].forward,
      face_metric_storage[1].inverse,
      face_metric_storage[1].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
    _fill_face_metric_storage_2d_axis2_kernel!(backend)(
      face_metric_storage[2].forward,
      face_metric_storage[2].inverse,
      face_metric_storage[2].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  elseif N == 3
    _fill_face_metric_storage_3d_axis1_kernel!(backend)(
      face_metric_storage[1].forward,
      face_metric_storage[1].inverse,
      face_metric_storage[1].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
    _fill_face_metric_storage_3d_axis2_kernel!(backend)(
      face_metric_storage[2].forward,
      face_metric_storage[2].inverse,
      face_metric_storage[2].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
    _fill_face_metric_storage_3d_axis3_kernel!(backend)(
      face_metric_storage[3].forward,
      face_metric_storage[3].inverse,
      face_metric_storage[3].conserved,
      edge,
      local_domain,
      global_domain,
      nhalo,
      t,
      params,
      T;
      ndrange=size(local_domain),
    )
  else
    throw(ArgumentError("Unsupported face metric dimension N=$N"))
  end
  KernelAbstractions.synchronize(backend)
  return nothing
end

#
# Legacy trait inference
#

function _coordinate_system_from_legacy(::OrthogonalGrid{D,T,CartesianCS}) where {D,T}
  CartesianCS()
end
function _coordinate_system_from_legacy(::OrthogonalGrid{D,T,CylindricalCS}) where {D,T}
  CylindricalCS()
end
function _coordinate_system_from_legacy(::OrthogonalGrid{D,T,SphericalCS}) where {D,T}
  SphericalCS()
end
_coordinate_system_from_legacy(::SphericalGrid1D) = SphericalCS()
_coordinate_system_from_legacy(::CylindricalGrid1D) = CylindricalCS()
_coordinate_system_from_legacy(::SphericalBasisCurvilinearGrid3D) = SphericalCS()

function _coordinate_system_from_legacy(mesh::AxisymmetricGrid2D)
  if mesh.rotational_axis === :x
    return AxisymmetricCS{:x}()
  else
    return AxisymmetricCS{:y}()
  end
end

function _coordinate_system_from_legacy(
  ::OrthogonalGrid{2,T,AxisymmetricCS{Axis}}
) where {T,Axis}
  Axis in (:x, :y) ||
    throw(ArgumentError("Unsupported axisymmetric axis `:$Axis` in legacy orthogonal grid"))
  AxisymmetricCS{Axis}()
end
_coordinate_system_from_legacy(::AbstractCurvilinearGrid) = CurvilinearCS()

_basis_trait_from_legacy(::SphericalBasisCurvilinearGrid3D) = SphericalBasis()
_basis_trait_from_legacy(::AbstractCurvilinearGrid) = CartesianBasis()
