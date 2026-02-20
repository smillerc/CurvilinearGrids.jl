#
# MappedGrid
#

struct MappedGrid{
  N,
  T,
  CS<:CoordinateSystemTrait,
  BT<:BasisTrait,
  NC,
  CC,
  MF,
  MFC,
  B,
  DB,
  DS,
  I,
  S,
  MC<:UnifiedMetricCaches,
} <: AbstractMappedOrDiscreteGrid
  node_coordinates::NC
  centroid_coordinates::CC
  mapping_functions::MF
  metric_functions_cache::MFC
  backend::B
  diff_backend::DB
  nhalo::Int
  discretization_scheme::DS
  discretization_scheme_name::Symbol
  iterators::I
  state::S
  metric_caches::MC
end

function _mapped_state(grid::MappedGrid)
  state = grid.state[]
  has_state = state isa NamedTuple && haskey(state, :t) && haskey(state, :params)
  has_state ? (state.t, state.params, true) : (nothing, nothing, false)
end

function _recompute_mapped_cell_metrics!(grid::MappedGrid{N,T}) where {N,T}
  t, params, has_state = _mapped_state(grid)
  if !has_state
    return nothing
  end
  cell_storage, _ = _ensure_metric_storage!(grid, Val(N), T)

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

function _recompute_mapped_face_metrics!(grid::MappedGrid{N,T}) where {N,T}
  t, params, has_state = _mapped_state(grid)
  if !has_state
    return nothing
  end
  _, face_storage = _ensure_metric_storage!(grid, Val(N), T)

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
  _recompute_mapped_cell_metrics!(grid)
  data = grid.metric_caches.cell.data

  if grid.metric_caches.cell.mode === :off
    return data
  end

  grid.metric_caches.cell.valid = true
  return data
end

function _refresh_face_metrics!(grid::MappedGrid; include_halo_region::Bool=false)
  _recompute_mapped_face_metrics!(grid)
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
  discretization_scheme::Symbol;
  backend,
  diff_backend,
  t,
  T::Type,
  compute_metrics::Bool,
  global_cell_indices,
  coordinate_system::CoordinateSystemTrait,
  basis::BasisTrait,
  cache_mode::Symbol,
) where {N}
  components = _build_unified_components(
    Val(N),
    mapping_functions,
    celldims,
    discretization_scheme,
    backend,
    diff_backend,
    T;
    global_cell_indices=global_cell_indices,
  )

  requested_mode = compute_metrics ? cache_mode : (cache_mode === :eager ? :lazy : cache_mode)
  caches = _new_metric_caches(
    requested_mode, components.cell_metric_storage, components.face_metric_storage
  )
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
    typeof(components.discretization_scheme),
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
    components.discretization_scheme,
    components.discretization_scheme_name,
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

  if requested_mode === :eager
    _refresh_cell_metrics!(grid)
    _refresh_face_metrics!(grid)
  end

  return grid
end

function MappedGrid(
  x::Function,
  params::NamedTuple,
  celldims::NTuple{1,Int},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  t=zero(Float64),
  T::Type=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  mapping_functions = (; x)
  return _new_mapped_grid(
    Val(1),
    mapping_functions,
    params,
    celldims,
    discretization_scheme;
    backend=backend,
    diff_backend=diff_backend,
    t=t,
    T=T,
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=cache_mode,
  )
end

function MappedGrid(
  x::Function,
  y::Function,
  params::NamedTuple,
  celldims::NTuple{2,Int},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  t=zero(Float64),
  T::Type=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  mapping_functions = (; x, y)
  return _new_mapped_grid(
    Val(2),
    mapping_functions,
    params,
    celldims,
    discretization_scheme;
    backend=backend,
    diff_backend=diff_backend,
    t=t,
    T=T,
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=cache_mode,
  )
end

function MappedGrid(
  x::Function,
  y::Function,
  z::Function,
  params::NamedTuple,
  celldims::NTuple{3,Int},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  t=zero(Float64),
  T::Type=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  mapping_functions = (; x, y, z)
  return _new_mapped_grid(
    Val(3),
    mapping_functions,
    params,
    celldims,
    discretization_scheme;
    backend=backend,
    diff_backend=diff_backend,
    t=t,
    T=T,
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=cache_mode,
  )
end

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
