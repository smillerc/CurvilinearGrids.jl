#
# MappedGrid
#

"""
    MappedGrid{N,T,CS,BT,...}

Unified grid backed by continuous mapping functions from computational space to
physical space.

`MappedGrid` stores node/centroid coordinates, mapping/metric function caches,
and independent metric caches for cell and face data.

# Fields
  - `node_coordinates`: Node coordinate arrays as `NTuple{N,AbstractArray}`.
  - `centroid_coordinates`: Cell-center coordinate arrays as `NTuple{N,AbstractArray}`.
  - `mapping_functions`: Physical mapping callbacks (`x1`, `x2`, `x3` by dimension).
  - `metric_functions_cache`: Cached metric-function closures.
  - `backend`: Compute backend used for storage allocation.
  - `diff_backend`: Differentiation backend used for metric evaluation.
  - `nhalo`: Halo width used by iterator/domain definitions.
  - `iterators`: Node/cell iterator bundle for local/global indexing.
  - `state`: Mutable reference to runtime state (`t`, `params`).
  - `metric_caches`: Independent cell and face metric caches.
"""
struct MappedGrid{N,T,CS<:CoordinateSystemTrait,BT<:BasisTrait,NC,CC,MF,MFC,B,DB,I,S,MC} <:
       AbstractMappedOrDiscreteGrid
  node_coordinates::NC
  centroid_coordinates::CC
  mapping_functions::MF
  metric_functions_cache::MFC
  backend::B
  diff_backend::DB
  nhalo::Int
  iterators::I
  state::S
  metric_caches::MC
end

function _mapped_state(grid::MappedGrid)
  state = grid.state[]
  has_state = state isa NamedTuple && haskey(state, :t) && haskey(state, :params)
  has_state ? (state.t, state.params, true) : (nothing, nothing, false)
end

function _recompute_mapped_cell_metrics!(
  grid::MappedGrid{N,T}; include_halo_region::Bool=false
) where {N,T}
  _require_metric_storage(grid, "refresh_cell_metrics!")
  t, params, has_state = _mapped_state(grid)
  if !has_state
    return nothing
  end
  cell_storage, _ = _ensure_metric_storage!(grid, Val(N), T)
  _fill_cell_metric_storage!(
    cell_storage,
    grid.metric_functions_cache,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
    T,
  )

  return nothing
end

function _recompute_mapped_face_metrics!(
  grid::MappedGrid{N,T}; include_halo_region::Bool=false
) where {N,T}
  _require_metric_storage(grid, "refresh_face_metrics!")
  t, params, has_state = _mapped_state(grid)
  if !has_state
    return nothing
  end
  _, face_storage = _ensure_metric_storage!(grid, Val(N), T)
  _fill_face_metric_storage!(
    face_storage,
    grid.metric_functions_cache,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
    T,
  )

  return nothing
end

function _refresh_cell_metrics!(grid::MappedGrid; include_halo_region::Bool=false)
  _require_metric_storage(grid, "refresh_cell_metrics!")
  _recompute_mapped_cell_metrics!(grid; include_halo_region=include_halo_region)
  data = grid.metric_caches.cell.data

  if grid.metric_caches.cell.mode === :off
    return data
  end

  grid.metric_caches.cell.valid = true
  return data
end

function _refresh_face_metrics!(grid::MappedGrid; include_halo_region::Bool=false)
  _require_metric_storage(grid, "refresh_face_metrics!")
  _recompute_mapped_face_metrics!(grid; include_halo_region=include_halo_region)
  data = grid.metric_caches.face.data

  if grid.metric_caches.face.mode === :off
    return data
  end

  grid.metric_caches.face.valid = true
  return data
end

function _new_mapped_grid(
  ::Val{N},
  mapping_functions,
  params::NamedTuple,
  celldims::NTuple{N,Int},
  nhalo::Int;
  backend,
  diff_backend,
  t,
  T::Type,
  compute_metrics::Bool,
  global_cell_indices,
  coordinate_system::CoordinateSystemTrait,
  basis::BasisTrait,
  cache_mode::Symbol,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait,
) where {N}
  _check_unified_basis_trait(basis)

  disable_metrics = (!compute_metrics && cache_mode === :off)
  components = _build_unified_components(
    Val(N),
    mapping_functions,
    celldims,
    nhalo,
    conserved_metric_scheme,
    backend,
    diff_backend,
    T;
    global_cell_indices=global_cell_indices,
    build_metric_storage=(!disable_metrics),
  )

  requested_mode =
    compute_metrics ? cache_mode : (cache_mode === :eager ? :lazy : cache_mode)
  caches = if disable_metrics
    nothing
  else
    _new_metric_caches(
      requested_mode, components.cell_metric_storage, components.face_metric_storage
    )
  end
  state = Ref((; t, params))

  grid = MappedGrid{
    N,
    T,
    typeof(coordinate_system),
    typeof(basis),
    typeof(components.node_coordinates),
    typeof(components.centroid_coordinates),
    typeof(mapping_functions),
    typeof(components.metric_functions_cache),
    typeof(backend),
    typeof(diff_backend),
    typeof(components.iterators),
    typeof(state),
    typeof(caches),
  }(
    components.node_coordinates,
    components.centroid_coordinates,
    mapping_functions,
    components.metric_functions_cache,
    backend,
    diff_backend,
    components.nhalo,
    components.iterators,
    state,
    caches,
  )

  _compute_unified_node_coordinates!(
    grid.node_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
  )
  _compute_unified_centroid_coordinates!(
    grid.centroid_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
  )

  if !disable_metrics && requested_mode === :eager
    _refresh_cell_metrics!(grid)
    _refresh_face_metrics!(grid)
  end

  return grid
end

"""
    MappedGrid(x[, y[, z]], params, celldims, nhalo; kwargs...)

Construct a mapped unified grid from continuous coordinate mapping functions.

# Arguments
  - `x`: Mapping function for the first physical coordinate.
  - `y`: Mapping function for the second physical coordinate (2D/3D).
  - `z`: Mapping function for the third physical coordinate (3D).
  - `params`: Mapping parameter tuple passed to mapping functions.
  - `celldims`: Cell counts in each computational dimension.
  - `nhalo`: Halo width used by node/cell domains.

# Keywords
  - `backend`: Storage backend. Default: `CPU()`.
  - `diff_backend`: Differentiation backend. Default: `AutoForwardDiff()`.
  - `t`: Initial time. Default: `zero(Float64)`.
  - `T`: Grid floating-point type. Default: `Float64`.
  - `compute_metrics`: Enable initial metric computation. Default: `true`.
    Set `compute_metrics=false` with `cache_mode=:off` to disable metric
    allocation/caching entirely.
  - `global_cell_indices`: Optional global index map. Default: `nothing`.
  - `coordinate_system`: Coordinate-system trait. Default: `CurvilinearCS()`.
  - `basis`: Basis trait. Default: `CartesianBasis()`.
  - `cache_mode`: Metric cache mode (`:eager`, `:lazy`, `:off`). Default: `:eager`.
  - `conserved_metric_scheme`: Conserved face interpolation scheme trait (`EdgeInterpolationOrder1()`, `EdgeInterpolationOrder2()`, `EdgeInterpolationOrder3()`). Default: `EdgeInterpolationOrder3()`.

# Returns
A `MappedGrid{N,T,...}` instance with initialized coordinates and metric cache
storage.
"""
function MappedGrid(
  x::Function,
  params::NamedTuple,
  celldims::NTuple{1,Int},
  nhalo::Integer;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  t=zero(Float64),
  T::Type=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  cache_mode::Symbol=:eager,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  nhalo_value = _normalize_nhalo(nhalo)
  mapping_functions = (; x1=x)
  return _new_mapped_grid(
    Val(1),
    mapping_functions,
    params,
    celldims,
    nhalo_value;
    backend=backend,
    diff_backend=diff_backend,
    t=t,
    T=T,
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=cache_mode,
    conserved_metric_scheme=conserved_metric_scheme,
  )
end

function MappedGrid(
  x::Function,
  y::Function,
  params::NamedTuple,
  celldims::NTuple{2,Int},
  nhalo::Integer;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  t=zero(Float64),
  T::Type=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  cache_mode::Symbol=:eager,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  nhalo_value = _normalize_nhalo(nhalo)
  mapping_functions = (; x1=x, x2=y)
  return _new_mapped_grid(
    Val(2),
    mapping_functions,
    params,
    celldims,
    nhalo_value;
    backend=backend,
    diff_backend=diff_backend,
    t=t,
    T=T,
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=cache_mode,
    conserved_metric_scheme=conserved_metric_scheme,
  )
end

function MappedGrid(
  x::Function,
  y::Function,
  z::Function,
  params::NamedTuple,
  celldims::NTuple{3,Int},
  nhalo::Integer;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  t=zero(Float64),
  T::Type=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  cache_mode::Symbol=:eager,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  nhalo_value = _normalize_nhalo(nhalo)
  mapping_functions = (; x1=x, x2=y, x3=z)
  return _new_mapped_grid(
    Val(3),
    mapping_functions,
    params,
    celldims,
    nhalo_value;
    backend=backend,
    diff_backend=diff_backend,
    t=t,
    T=T,
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=cache_mode,
    conserved_metric_scheme=conserved_metric_scheme,
  )
end

"""
    update!(grid::MappedGrid{N}, t::Real, params::NamedTuple) where {N}

Update mapped grid coordinates and invalidate metric caches for a new state.

# Arguments
  - `grid`: Target mapped grid.
  - `t`: New time value.
  - `params`: New parameter tuple used by mapping functions.

# Returns
`nothing`.
"""
function update!(grid::MappedGrid{N}, t::Real, params::NamedTuple) where {N}
  _compute_unified_node_coordinates!(
    grid.node_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
  )
  _compute_unified_centroid_coordinates!(
    grid.centroid_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
  )

  grid.state[] = (; t, params)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end

function update!(grid::MappedGrid, t::Real=zero(Float64))
  state = grid.state[]
  params = state isa NamedTuple && haskey(state, :params) ? state.params : (;)
  return update!(grid, t, params)
end
