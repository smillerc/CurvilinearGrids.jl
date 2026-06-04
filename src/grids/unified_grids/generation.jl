# Unified grid coordinate, storage, and volume generation.

function get_iterators(celldims::NTuple{N,Int}, nhalo::Int, global_cell_domain) where {N}
  cellCI = CartesianIndices(celldims .+ 2nhalo)
  nodeCI = CartesianIndices(celldims .+ 1 .+ 2nhalo)

  node = (full=nodeCI, domain=expand(nodeCI, -nhalo))
  cell = (full=cellCI, domain=expand(cellCI, -nhalo))

  if isnothing(global_cell_domain)
    global_domain = (node=node, cell=cell)
  else
    global_domain = (node=expand_upper(global_cell_domain, +1), cell=global_cell_domain)
  end
  return (; node, cell, nhalo, global_domain)
end


# Metric and coordinate storage allocation.

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


# Coordinate storage allocation.

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



# Unified component construction.

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


# Coordinate generation kernels and launchers.

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

function _compute_unified_face_coordinates!(
  face_coordinates, mapping_functions, iterators, nhalo::Int, t, params, ::Val{N}, backend
) where {N}
  _compute_unified_face_coordinates!(
    backend, face_coordinates, mapping_functions, iterators, nhalo, t, params, Val(N)
  )
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

# Cell geometry and per-cell volume helpers.

"""
    cellvolume(grid, idx)

Compute cell volume at a given index.

For mapped/discrete grids, `idx` may be a `CartesianIndex`, an integer tuple for
stored cell metrics, or a real-valued tuple for continuous-coordinate metric
evaluation. For orthogonal grids, `idx` is an integer cell index.

# Arguments
  - `grid`: Unified grid instance.
  - `idx`: `CartesianIndex`, integer tuple, or supported real tuple.
"""
function cellvolume(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex)
  cellvolume(grid, idx.I)
end
function cellvolume(
  grid::Union{MappedGrid{N,T,CS,BT},DiscreteGrid{N,T,CS,BT}}, idx::NTuple{N,Int}
) where {N,T,CS,BT}
  _cellvolume_dispatch(CS(), BT(), grid, idx)
end
function cellvolume(
  grid::Union{MappedGrid{N,T,CS,BT},DiscreteGrid{N,T,CS,BT}}, idx::Tuple{Vararg{Real,N}}
) where {N,T,CS,BT}
  _cellvolume_dispatch(CS(), BT(), grid, _promote_real_tuple(idx))
end

"""
    cell_jacobian(grid, idx)

Return the geometric Jacobian factor used by diffusion operators at cell index
`idx`.

For mapped/discrete grids this is `abs(forward_cell_metrics(grid, idx).J)`.
For orthogonal grids it is the reduced volumetric measure associated with the
coordinate family, for example `r`, `r^2`, or `r^2 sin(θ)`.
"""
@inline function cell_jacobian(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex)
  cell_jacobian(grid, idx.I)
end
@inline function cell_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  abs(forward_cell_metrics(grid, idx).J)
end

@inline function _orth_cell_coord(
  grid::OrthogonalGrid{N}, dim::Int, idx::NTuple{N,Int}
) where {N}
  grid.centroid_coordinates[dim][idx[dim]]
end
@inline function _orth_face_coord(
  grid::OrthogonalGrid{N}, dim::Int, idx::NTuple{N,Int}
) where {N}
  grid.node_coordinates[dim][idx[dim]]
end

@inline function _axisymmetric_radial_dim(::AxisymmetricCS{:y}, ::Val{2})
  1
end
@inline function _axisymmetric_radial_dim(::AxisymmetricCS{:x}, ::Val{2})
  2
end
function _axisymmetric_radial_dim(::AxisymmetricCS{Axis}, ::Val{2}) where {Axis}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis`. Supported values are `:x` and `:y`."
    ),
  )
end

@inline function cell_jacobian(grid::AbstractOrthogonalGrid, idx::CartesianIndex)
  cell_jacobian(grid, idx.I)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{N,T,CartesianCS}, idx::NTuple{N,Int}
) where {N,T}
  one(T)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{1,T,CylindricalCS}, idx::NTuple{1,Int}
) where {T}
  _orth_cell_coord(grid, 1, idx)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{2,T,AxisymmetricCS{Axis}}, idx::NTuple{2,Int}
) where {T,Axis}
  ridx = _axisymmetric_radial_dim(AxisymmetricCS{Axis}(), Val(2))
  _orth_cell_coord(grid, ridx, idx)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{1,T,SphericalCS}, idx::NTuple{1,Int}
) where {T}
  r = _orth_cell_coord(grid, 1, idx)
  r^2
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{2,T,SphericalCS}, idx::NTuple{2,Int}
) where {T}
  r = _orth_cell_coord(grid, 1, idx)
  θ = _orth_cell_coord(grid, 2, idx)
  r^2 * sin(θ)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{3,T,SphericalCS}, idx::NTuple{3,Int}
) where {T}
  r = _orth_cell_coord(grid, 1, idx)
  θ = _orth_cell_coord(grid, 2, idx)
  r^2 * sin(θ)
end
function cell_jacobian(grid::AbstractOrthogonalGrid, idx::NTuple{N,Int}) where {N}
  throw(
    ArgumentError(
      "`cell_jacobian` is undefined for $(typeof(grid)) with $(N)-D integer indices."
    ),
  )
end


# Cell volume array generation and grid sizes.

@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  return abs(_cell_forward_metric_at(grid, idx).J)
end
@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  return abs(T(det(_continuous_forward_jacobian(grid, idx))))
end

@inline _radial_centroid_1d(
  grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::NTuple{1,Int}
) = grid.centroid_coordinates[1][idx...]
@inline _radial_centroid_1d(
  grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::Tuple{Vararg{Real,1}}
) = _continuous_coord(grid, idx)[1]

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:x}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::NTuple{2,Int}
)
  return grid.centroid_coordinates[2][idx...]
end
@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:x},
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  return _continuous_coord(grid, idx)[2]
end

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:y}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::NTuple{2,Int}
)
  return grid.centroid_coordinates[1][idx...]
end
@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:y},
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  return _continuous_coord(grid, idx)[1]
end

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{Axis}, grid::AbstractMappedOrDiscreteGrid, idx::Tuple{Vararg{Real,2}}
) where {Axis}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis`. Supported values are `:x` and `:y`."
    ),
  )
end

@inline function _cellvolume_dispatch(
  ::CoordinateSystemTrait,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::NTuple{N,Int},
) where {N}
  return _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CoordinateSystemTrait,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,N}},
) where {N}
  return _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS, ::CartesianBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{1,Int}
)
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::NTuple{2,Int},
)
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::NTuple{2,Int},
) where {Axis}
  return conservation_cell_metric_scale(cs, CartesianBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
) where {Axis}
  return conservation_cell_metric_scale(cs, CartesianBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{1,Int}
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::NTuple{3,Int},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::Tuple{Vararg{Real,3}},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{2,Int}
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  cs::CoordinateSystemTrait,
  bt::SphericalBasis,
  ::AbstractMappedOrDiscreteGrid,
  ::NTuple{N,Int},
) where {N}
  throw(
    ArgumentError(
      "Unsupported basis/coordinate-system combination: $(typeof(bt)) with $(typeof(cs))."
    ),
  )
end

function _cellvolume_dispatch(
  cs::CoordinateSystemTrait,
  bt::SphericalBasis,
  ::AbstractMappedOrDiscreteGrid,
  ::Tuple{Vararg{Real,N}},
) where {N}
  throw(
    ArgumentError(
      "Unsupported basis/coordinate-system combination: $(typeof(bt)) with $(typeof(cs))."
    ),
  )
end

function cellvolumes(grid::Union{MappedGrid{N},DiscreteGrid{N}}; include_halo::Bool=false) where {N}
  volumes = KernelAbstractions.zeros(
    grid.backend, eltype(grid), size(grid.iterators.cell.full); unified=false
  )
  _compute_unified_cell_volumes!(volumes, grid, Val(N))
  include_halo && return volumes
  @views return volumes[grid.iterators.cell.domain]
end

function cellvolumes(grid::OrthogonalGrid; include_halo::Bool=false)
  include_halo && return grid.cell_volumes
  @views return grid.cell_volumes[grid.iterators.cell.domain]
end

function _compute_unified_cell_volumes!(volumes, grid, ::Val{N}) where {N}
  backend = KernelAbstractions.get_backend(volumes)
  kernel = _compute_unified_cell_volumes_kernel!(backend)
  metrics = cell_metrics(grid)
  kernel(
    volumes,
    metrics.forward,
    grid.centroid_coordinates,
    coordinate_system(grid),
    basis_trait(grid),
    grid.iterators.cell.full,
    Val(N);
    ndrange=length(grid.iterators.cell.full),
  )
  KernelAbstractions.synchronize(backend)
  return volumes
end

@kernel function _compute_unified_cell_volumes_kernel!(
  volumes,
  cell_forward_metrics,
  centroid_coordinates,
  coordinate_system,
  basis,
  local_domain,
  ::Val{N},
) where {N}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  q = SVector{N}(ntuple(d -> centroid_coordinates[d][I], N))
  scale = conservation_cell_metric_scale(coordinate_system, basis, q)
  volumes[I] = scale * abs(cell_forward_metrics[I].J)
end

cellsize(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.domain)
cellsize(grid::OrthogonalGrid) = size(grid.iterators.cell.domain)

cellsize_withhalo(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.full)
cellsize_withhalo(grid::OrthogonalGrid) = size(grid.iterators.cell.full)
