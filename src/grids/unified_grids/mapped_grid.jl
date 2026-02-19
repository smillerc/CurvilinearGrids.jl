#
# MappedGrid
#

mutable struct MappedGrid{C,CS<:CoordinateSystemTrait,BT<:BasisTrait} <:
               AbstractMappedOrDiscreteGrid
  core::C
  coordinate_system_trait::CS
  basis_vector_trait::BT
  state::Any
  metric_caches::UnifiedMetricCaches
end

function _mapped_state(grid::MappedGrid)
  state = grid.state
  has_state = state isa NamedTuple && haskey(state, :t) && haskey(state, :params)

  if has_state
    return state.t, state.params, true
  else
    return nothing, nothing, false
  end
end

function _recompute_mapped_cell_metrics!(grid::MappedGrid)
  t, params, has_state = _mapped_state(grid)
  if !has_state
    return nothing
  end

  # Pull in the continuous-grid metric pipeline directly.
  compute_node_coordinates!(grid.core, t, params)
  compute_centroid_coordinates!(grid.core, t, params)
  compute_cell_metrics!(grid.core, t, params)
  return nothing
end

function _recompute_mapped_face_metrics!(grid::MappedGrid)
  t, params, has_state = _mapped_state(grid)
  if !has_state
    return nothing
  end

  # Face metrics depend on the current cell metrics.
  _recompute_mapped_cell_metrics!(grid)
  compute_edge_metrics!(grid.core, t, params)
  return nothing
end

function _refresh_cell_metrics!(grid::MappedGrid; include_halo_region::Bool=false)
  _recompute_mapped_cell_metrics!(grid)
  data = getproperty(grid.core, :cell_center_metrics)

  if grid.metric_caches.cell.mode === :off
    return data
  end

  grid.metric_caches.cell.data = deepcopy(data)
  grid.metric_caches.cell.valid = true
  return grid.metric_caches.cell.data
end

function _refresh_face_metrics!(grid::MappedGrid; include_halo_region::Bool=false)
  _recompute_mapped_face_metrics!(grid)
  data = getproperty(grid.core, :edge_metrics)

  if grid.metric_caches.face.mode === :off
    return data
  end

  grid.metric_caches.face.data = deepcopy(data)
  grid.metric_caches.face.valid = true
  return grid.metric_caches.face.data
end

function MappedGrid(
  core::Union{
    AbstractContinuousCurvilinearGrid1D,
    AbstractContinuousCurvilinearGrid2D,
    AbstractContinuousCurvilinearGrid3D,
  };
  coordinate_system::CoordinateSystemTrait=_coordinate_system_from_legacy(core),
  basis::BasisTrait=_basis_trait_from_legacy(core),
  state=nothing,
  cache_mode::Symbol=:eager,
)
  caches = _new_metric_caches(cache_mode)
  grid = MappedGrid(core, coordinate_system, basis, state, caches)

  if cache_mode === :eager
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
  T=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  core = ContinuousCurvilinearGrid1D(
    x,
    params,
    celldims,
    discretization_scheme,
    backend,
    diff_backend,
    t,
    T;
    compute_metrics=false,
    global_cell_indices=global_cell_indices,
  )

  requested_mode = compute_metrics ? cache_mode : (cache_mode === :eager ? :lazy : cache_mode)
  return MappedGrid(
    core;
    coordinate_system=coordinate_system,
    basis=basis,
    state=(; t, params),
    cache_mode=requested_mode,
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
  T=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  core = ContinuousCurvilinearGrid2D(
    x,
    y,
    params,
    celldims,
    discretization_scheme,
    backend,
    diff_backend,
    t,
    T;
    compute_metrics=false,
    global_cell_indices=global_cell_indices,
  )

  requested_mode = compute_metrics ? cache_mode : (cache_mode === :eager ? :lazy : cache_mode)
  return MappedGrid(
    core;
    coordinate_system=coordinate_system,
    basis=basis,
    state=(; t, params),
    cache_mode=requested_mode,
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
  T=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  core = ContinuousCurvilinearGrid3D(
    x,
    y,
    z,
    params,
    celldims,
    discretization_scheme,
    backend,
    diff_backend,
    t,
    T;
    compute_metrics=false,
    global_cell_indices=global_cell_indices,
  )

  requested_mode = compute_metrics ? cache_mode : (cache_mode === :eager ? :lazy : cache_mode)
  return MappedGrid(
    core;
    coordinate_system=coordinate_system,
    basis=basis,
    state=(; t, params),
    cache_mode=requested_mode,
  )
end

function update!(grid::MappedGrid, t, params)
  # Keep coordinate fields in sync immediately; metric caches refresh independently.
  compute_node_coordinates!(grid.core, t, params)
  compute_centroid_coordinates!(grid.core, t, params)

  grid.state = (; t, params)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end

function update!(grid::MappedGrid, args...; kwargs...)
  update!(grid.core, args...; kwargs...)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end
