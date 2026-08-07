# MappedGrid type and constructors.

struct _ReparameterizedMappingComponent{D,F,O,S} <: Function
  f::F
  origin::O
  scale::S
end

@inline _reparameterized_coordinate(q, origin, scale) = origin + (q - oneunit(q)) * scale

@inline _mapping_component_keys(::Val{1}) = (:x1,)
@inline _mapping_component_keys(::Val{2}) = (:x1, :x2)
@inline _mapping_component_keys(::Val{3}) = (:x1, :x2, :x3)

@inline function (m::_ReparameterizedMappingComponent{1})(t, ξ, params)
  return m.f(t, _reparameterized_coordinate(ξ, m.origin[1], m.scale[1]), params)
end

@inline function (m::_ReparameterizedMappingComponent{2})(t, ξ, η, params)
  ξf = _reparameterized_coordinate(ξ, m.origin[1], m.scale[1])
  ηf = _reparameterized_coordinate(η, m.origin[2], m.scale[2])
  return m.f(t, ξf, ηf, params)
end

@inline function (m::_ReparameterizedMappingComponent{3})(t, ξ, η, ζ, params)
  ξf = _reparameterized_coordinate(ξ, m.origin[1], m.scale[1])
  ηf = _reparameterized_coordinate(η, m.origin[2], m.scale[2])
  ζf = _reparameterized_coordinate(ζ, m.origin[3], m.scale[3])
  return m.f(t, ξf, ηf, ζf, params)
end

"""
    reparameterize_mapping(mapping_functions, origin, scale)

Return mapping components that evaluate `mapping_functions` after the affine
one-based computational-coordinate transform
`q_global[d] = origin[d] + (q_local[d] - 1) * scale[d]`.

`mapping_functions` must be a named tuple with exactly the keys `:x1`, `:x2`,
and `:x3` required by its dimension. The returned named tuple preserves those
keys and contains concrete callable values suitable for mapped-grid CPU or GPU
evaluation. `origin` and `scale` are stored by value; no mapping parameters or
mutable grid state are captured.
"""
function reparameterize_mapping(
  mapping_functions::NamedTuple, origin::NTuple{D,<:Real}, scale::NTuple{D,<:Real}
) where {D}
  D in 1:3 || throw(
    ArgumentError("Mapping reparameterization supports dimensions 1, 2, and 3; got D=$D.")
  )
  expected_keys = _mapping_component_keys(Val(D))
  keys(mapping_functions) == expected_keys || throw(
    ArgumentError(
      "Expected mapping-function keys $expected_keys for D=$D; got $(keys(mapping_functions)).",
    ),
  )
  components = ntuple(Val(D)) do d
    f = getfield(mapping_functions, d)
    _ReparameterizedMappingComponent{D,typeof(f),typeof(origin),typeof(scale)}(
      f, origin, scale
    )
  end
  return NamedTuple{expected_keys}(components)
end

function reparameterize_mapping(
  mapping_functions::NamedTuple, origin::Tuple{Vararg{Real}}, scale::Tuple{Vararg{Real}}
)
  throw(
    DimensionMismatch(
      "Mapping origin and scale dimensions must match; got $(length(origin)) and $(length(scale)).",
    ),
  )
end

"""
    MappedGrid{N,T,CS,BT,...}

Unified grid backed by continuous mapping functions from computational space to
physical space.

`MappedGrid` stores node/centroid/face coordinates, mapping/metric function caches,
and independent metric caches for cell and face data.

# Fields
  - `node_coordinates`: Node coordinate arrays as `NTuple{N,AbstractArray}`.
  - `centroid_coordinates`: Cell-center coordinate arrays as `NTuple{N,AbstractArray}`.
  - `face_coordinates`: High-side face-center coordinate arrays indexed by face axis.
  - `mapping_functions`: Physical mapping callbacks (`x1`, `x2`, `x3` by dimension).
  - `metric_functions_cache`: Cached metric-function closures.
  - `backend`: Compute backend used for storage allocation.
  - `diff_backend`: Differentiation backend used for metric evaluation.
  - `nhalo`: Halo width used by iterator/domain definitions.
  - `iterators`: Node/cell iterator bundle for local/global indexing.
  - `state`: Mutable reference to runtime state (`t`, `params`).
  - `metric_caches`: Independent cell and face metric caches.
  - `conserved_metric_scheme`: Scheme used to construct conserved face metrics.
"""
struct MappedGrid{
  N,
  T,
  CS<:CoordinateSystemTrait,
  BT<:BasisTrait,
  NC,
  CC,
  FC,
  MF,
  MFC,
  B,
  DB,
  I,
  S,
  MC,
  CMS<:EdgeInterpolationSchemeTrait,
} <: AbstractMappedOrDiscreteGrid
  node_coordinates::NC
  centroid_coordinates::CC
  face_coordinates::FC
  mapping_functions::MF
  metric_functions_cache::MFC
  backend::B
  diff_backend::DB
  nhalo::Int
  iterators::I
  state::S
  metric_caches::MC
  conserved_metric_scheme::CMS
end

function Base.show(io::IO, grid_type::Type{<:MappedGrid{N,T}}) where {N,T}
  if get(io, :backtrace, false)::Bool
    print(io, "MappedGrid{", N, ", ", T, ", …}")
    return nothing
  end
  return invoke(Base.show, Tuple{IO,Type}, io, grid_type)
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
    grid.backend,
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
    grid.backend,
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
  metric_functions_cache=nothing,
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
    metric_functions_cache=metric_functions_cache,
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
    typeof(components.face_coordinates),
    typeof(mapping_functions),
    typeof(components.metric_functions_cache),
    typeof(backend),
    typeof(diff_backend),
    typeof(components.iterators),
    typeof(state),
    typeof(caches),
    typeof(conserved_metric_scheme),
  }(
    components.node_coordinates,
    components.centroid_coordinates,
    components.face_coordinates,
    mapping_functions,
    components.metric_functions_cache,
    backend,
    diff_backend,
    components.nhalo,
    components.iterators,
    state,
    caches,
    conserved_metric_scheme,
  )

  _compute_unified_node_coordinates!(
    grid.node_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
    grid.backend,
  )
  _compute_unified_centroid_coordinates!(
    grid.centroid_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
    grid.backend,
  )
  _compute_unified_face_coordinates!(
    grid.face_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
    grid.backend,
  )

  if !disable_metrics && requested_mode === :eager
    _refresh_cell_metrics!(grid)
    _refresh_face_metrics!(grid)
  end

  return grid
end

"""
    local_grid(grid::MappedGrid, global_cell_domain; kwargs...)

Materialize a rank-local `MappedGrid` over the supplied block-global, non-halo
cell domain. The analytic mapping and metric functions are preserved, while
coordinate and metric storage are allocated only for the local cells and halos.

`global_cell_domain` may be a `CartesianIndices` domain or a tuple of integer
ranges. Metrics are computed after localization by default.
"""
function local_grid(
  grid::MappedGrid{N,T},
  global_cell_domain::CartesianIndices{N};
  backend=grid.backend,
  compute_metrics::Bool=true,
  cache_mode::Symbol=:eager,
) where {N,T}
  state = grid.state[]
  has_state = state isa NamedTuple && haskey(state, :t) && haskey(state, :params)
  has_state || throw(ArgumentError("MappedGrid has no valid mapping state to localize."))

  return _new_mapped_grid(
    Val(N),
    grid.mapping_functions,
    state.params,
    size(global_cell_domain),
    grid.nhalo;
    backend=backend,
    diff_backend=grid.diff_backend,
    t=state.t,
    T=T,
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_domain,
    coordinate_system=coordinate_system(grid),
    basis=basis_trait(grid),
    cache_mode=cache_mode,
    conserved_metric_scheme=CurvatureCorrectedReconstruction(),
    metric_functions_cache=grid.metric_functions_cache,
  )
end

function local_grid(
  grid::MappedGrid{N}, global_cell_ranges::NTuple{N,<:AbstractUnitRange}; kwargs...
) where {N}
  return local_grid(grid, CartesianIndices(global_cell_ranges); kwargs...)
end

"""
    MappedGrid(x1[, x2[, x3]], params, celldims, nhalo; kwargs...)

Construct a mapped unified grid from continuous coordinate mapping functions.

Mapping callbacks are evaluated as `x1(t, xi, params)` in 1D,
`x1(t, xi, eta, params)`/`x2(t, xi, eta, params)` in 2D, and with
`zeta` added in 3D.

# Arguments
  - `x1`: Mapping function for the first physical coordinate.
  - `x2`: Mapping function for the second physical coordinate (2D/3D).
  - `x3`: Mapping function for the third physical coordinate (3D).
  - `params`: Mapping parameter tuple passed to mapping functions.
  - `celldims`: Cell counts in each computational dimension.
  - `nhalo`: Halo width used by node/cell domains.

# Keywords
  - `backend`: Storage backend. Default: `CPU()`.
  - `diff_backend`: Differentiation backend. Default: `AutoForwardDiff()`.
  - `t`: Initial time. Default: `zero(Float64)`.
  - `T`: Grid floating-point type. Default: `Float64`.
  - `compute_metrics`: Enable initial metric-cache computation. Default: `true`.
    When `false`, metrics are computed on first access unless `cache_mode=:off`;
    set `compute_metrics=false, cache_mode=:off` to disable metric storage
    entirely.
  - `global_cell_indices`: Optional block-global `CartesianIndices` domain for
    the non-halo cells. Its size must equal `celldims`. Default: `nothing`.
  - `coordinate_system`: Coordinate-system trait. Default: `CurvilinearCS()`.
  - `basis`: Basis trait. Default: `CartesianBasis()`.
  - `cache_mode`: Metric cache mode (`:eager`, `:lazy`, `:off`). Default: `:eager`.
  - `conserved_metric_scheme`: Conserved face reconstruction selector
    (`EndpointAverageReconstruction()`, `GradientCorrectedReconstruction()`,
    `CurvatureCorrectedReconstruction()`). Default:
    `CurvatureCorrectedReconstruction()`.

# Returns
A `MappedGrid{N,T,...}` instance with initialized coordinates and metric cache
storage.
"""
function MappedGrid(
  x1::Function,
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
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  @assert nhalo >= 0

  mapping_functions = (; x1=x1)
  return _new_mapped_grid(
    Val(1),
    mapping_functions,
    params,
    celldims,
    nhalo;
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
  x1::Function,
  x2::Function,
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
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  @assert nhalo >= 0

  mapping_functions = (; x1=x1, x2=x2)
  return _new_mapped_grid(
    Val(2),
    mapping_functions,
    params,
    celldims,
    nhalo;
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
  x1::Function,
  x2::Function,
  x3::Function,
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
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=CurvatureCorrectedReconstruction(),
)
  @assert nhalo >= 0

  mapping_functions = (; x1=x1, x2=x2, x3=x3)
  return _new_mapped_grid(
    Val(3),
    mapping_functions,
    params,
    celldims,
    nhalo;
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
    grid.backend,
  )
  _compute_unified_centroid_coordinates!(
    grid.centroid_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
    grid.backend,
  )
  _compute_unified_face_coordinates!(
    grid.face_coordinates,
    grid.mapping_functions,
    grid.iterators,
    grid.nhalo,
    t,
    params,
    Val(N),
    grid.backend,
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
