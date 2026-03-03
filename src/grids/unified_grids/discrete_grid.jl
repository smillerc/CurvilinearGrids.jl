#
# DiscreteGrid
#

using Interpolations

"""
    DiscreteGrid{N,T,CS,BT,...}

Unified grid built from user-provided coordinate arrays and linear
interpolation in computational space.

`DiscreteGrid` stores interpolation-backed mapping functions along with
node/centroid coordinates and independent cell/face metric caches.

# Fields
  - `node_coordinates`: Node coordinate arrays as `NTuple{N,AbstractArray}`.
  - `centroid_coordinates`: Cell-center coordinate arrays as `NTuple{N,AbstractArray}`.
  - `mapping_functions`: Interpolation-backed mapping callbacks.
  - `metric_functions_cache`: Cached metric-function closures.
  - `backend`: Compute backend used for storage allocation.
  - `diff_backend`: Differentiation backend used for metric evaluation.
  - `nhalo`: Halo width used by iterator/domain definitions.
  - `iterators`: Node/cell iterator bundle for local/global indexing.
  - `interpolation`: Interpolation mode symbol (`:linear`).
  - `interpolants`: Cached interpolation objects.
  - `state`: Mutable reference to runtime state (`t`, `params`).
  - `metric_caches`: Independent cell and face metric caches.
"""
struct DiscreteGrid{
  N,T,CS<:CoordinateSystemTrait,BT<:BasisTrait,IP,NC,CC,MF,MFC,B,DB,I,S,MC
} <: AbstractMappedOrDiscreteGrid
  node_coordinates::NC
  centroid_coordinates::CC
  mapping_functions::MF
  metric_functions_cache::MFC
  backend::B
  diff_backend::DB
  nhalo::Int
  iterators::I
  interpolation::Symbol
  interpolants::IP
  state::S
  metric_caches::MC
end

function _strip_halo_nodes(A, nhalo::Int)
  if nhalo <= 0
    return A
  end
  inds = ntuple(d -> (firstindex(A, d) + nhalo):(lastindex(A, d) - nhalo), ndims(A))
  return @view A[inds...]
end

_linear_interpolant(A) = extrapolate(interpolate(A, BSpline(Linear())), Line())

struct InterpolantMapping{N,I}
  itp::I
end

@inline (m::InterpolantMapping{1})(_t, ξ, _p) = m.itp(ξ)
@inline (m::InterpolantMapping{2})(_t, ξ, η, _p) = m.itp(ξ, η)
@inline (m::InterpolantMapping{3})(_t, ξ, η, ζ, _p) = m.itp(ξ, η, ζ)

function MetricCache(
  x::InterpolantMapping{1},
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  return MetricCache(
    (t, ξ, p) -> x(t, ξ, p), backend; edge_interpolation_scheme=edge_interpolation_scheme
  )
end

function MetricCache(
  x::InterpolantMapping{2},
  y::InterpolantMapping{2},
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  return MetricCache(
    (t, ξ, η, p) -> x(t, ξ, η, p),
    (t, ξ, η, p) -> y(t, ξ, η, p),
    backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )
end

function MetricCache(
  x::InterpolantMapping{3},
  y::InterpolantMapping{3},
  z::InterpolantMapping{3},
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  return MetricCache(
    (t, ξ, η, ζ, p) -> x(t, ξ, η, ζ, p),
    (t, ξ, η, ζ, p) -> y(t, ξ, η, ζ, p),
    (t, ξ, η, ζ, p) -> z(t, ξ, η, ζ, p),
    backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )
end

@inline function _discrete_jacobian_from_interpolants(interpolants::NamedTuple, ::Val{1})
  function jacobian(_t, ξ, _p)
    gx = Interpolations.gradient(interpolants.x1, ξ)
    return @SMatrix [gx[1]]
  end
  return jacobian
end

@inline function _discrete_jacobian_from_interpolants(interpolants::NamedTuple, ::Val{2})
  function jacobian(_t, ξ, η, _p)
    gx = Interpolations.gradient(interpolants.x1, ξ, η)
    gy = Interpolations.gradient(interpolants.x2, ξ, η)
    return @SMatrix [gx[1] gx[2]; gy[1] gy[2]]
  end
  return jacobian
end

@inline function _discrete_jacobian_from_interpolants(interpolants::NamedTuple, ::Val{3})
  function jacobian(_t, ξ, η, ζ, _p)
    gx = Interpolations.gradient(interpolants.x1, ξ, η, ζ)
    gy = Interpolations.gradient(interpolants.x2, ξ, η, ζ)
    gz = Interpolations.gradient(interpolants.x3, ξ, η, ζ)
    return @SMatrix [
      gx[1] gx[2] gx[3]
      gy[1] gy[2] gy[3]
      gz[1] gz[2] gz[3]
    ]
  end
  return jacobian
end

@inline function _discrete_metric_cache(
  ::Val{1}, metric_cache::MetricCache, interpolants::NamedTuple
)
  jacobian = _discrete_jacobian_from_interpolants(interpolants, Val(1))
  J(t, ξ, p) = det(jacobian(t, ξ, p))
  forward = merge(metric_cache.forward, (; jacobian=jacobian, J=J))
  return MetricCache(forward, metric_cache.inverse, metric_cache.edge)
end

@inline function _discrete_metric_cache(
  ::Val{2}, metric_cache::MetricCache, interpolants::NamedTuple
)
  jacobian = _discrete_jacobian_from_interpolants(interpolants, Val(2))
  J(t, ξ, η, p) = det(jacobian(t, ξ, η, p))
  forward = merge(metric_cache.forward, (; jacobian=jacobian, J=J))
  return MetricCache(forward, metric_cache.inverse, metric_cache.edge)
end

@inline function _discrete_metric_cache(
  ::Val{3}, metric_cache::MetricCache, interpolants::NamedTuple
)
  jacobian = _discrete_jacobian_from_interpolants(interpolants, Val(3))
  J(t, ξ, η, ζ, p) = det(jacobian(t, ξ, η, ζ, p))
  forward = merge(metric_cache.forward, (; jacobian=jacobian, J=J))
  return MetricCache(forward, metric_cache.inverse, metric_cache.edge)
end

function _validate_discrete_interpolation(interpolation::Symbol)
  if interpolation !== :linear
    throw(
      ArgumentError(
        "`DiscreteGrid` currently supports only linear interpolation (`interpolation=:linear`).",
      ),
    )
  end
  return interpolation
end

function _discrete_state(grid::DiscreteGrid)
  state = grid.state[]
  has_state = state isa NamedTuple && haskey(state, :t) && haskey(state, :params)
  has_state ? (state.t, state.params, true) : (nothing, nothing, false)
end

function _recompute_discrete_cell_metrics!(
  grid::DiscreteGrid{N,T}; include_halo_region::Bool=false
) where {N,T}
  _require_metric_storage(grid, "refresh_cell_metrics!")
  t, params, has_state = _discrete_state(grid)
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

function _recompute_discrete_face_metrics!(
  grid::DiscreteGrid{N,T}; include_halo_region::Bool=false
) where {N,T}
  _require_metric_storage(grid, "refresh_face_metrics!")
  t, params, has_state = _discrete_state(grid)
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

function _refresh_cell_metrics!(grid::DiscreteGrid; include_halo_region::Bool=false)
  _require_metric_storage(grid, "refresh_cell_metrics!")
  _recompute_discrete_cell_metrics!(grid; include_halo_region=include_halo_region)
  data = grid.metric_caches.cell.data

  if grid.metric_caches.cell.mode === :off
    return data
  end

  grid.metric_caches.cell.valid = true
  return data
end

function _refresh_face_metrics!(grid::DiscreteGrid; include_halo_region::Bool=false)
  _require_metric_storage(grid, "refresh_face_metrics!")
  _recompute_discrete_face_metrics!(grid; include_halo_region=include_halo_region)
  data = grid.metric_caches.face.data

  if grid.metric_caches.face.mode === :off
    return data
  end

  grid.metric_caches.face.valid = true
  return data
end

function _new_discrete_grid(
  ::Val{N},
  mapping_functions,
  interpolants,
  celldims::NTuple{N,Int},
  nhalo::Int;
  backend,
  diff_backend,
  t,
  params,
  T::Type,
  compute_metrics::Bool,
  coordinate_system::CoordinateSystemTrait,
  basis::BasisTrait,
  interpolation::Symbol,
  cache_mode::Symbol,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait,
) where {N}
  _check_unified_basis_trait(basis)
  _validate_discrete_interpolation(interpolation)

  @assert nhalo >= 0

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
    build_metric_storage=(!disable_metrics),
  )
  discrete_metric_functions_cache = _discrete_metric_cache(
    Val(N), components.metric_functions_cache, interpolants
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
  grid = DiscreteGrid{
    N,
    T,
    typeof(coordinate_system),
    typeof(basis),
    typeof(interpolants),
    typeof(components.node_coordinates),
    typeof(components.centroid_coordinates),
    typeof(mapping_functions),
    typeof(discrete_metric_functions_cache),
    typeof(backend),
    typeof(diff_backend),
    typeof(components.iterators),
    typeof(state),
    typeof(caches),
  }(
    components.node_coordinates,
    components.centroid_coordinates,
    mapping_functions,
    discrete_metric_functions_cache,
    backend,
    diff_backend,
    components.nhalo,
    components.iterators,
    interpolation,
    interpolants,
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

  if !disable_metrics && requested_mode === :eager
    _refresh_cell_metrics!(grid)
    _refresh_face_metrics!(grid)
  end

  return grid
end

"""
    DiscreteGrid(x[, y[, z]], nhalo; kwargs...)

Construct a discrete unified grid from coordinate arrays using linear
interpolation.

# Arguments
  - `x`: First coordinate array.
  - `y`: Second coordinate array (2D/3D).
  - `z`: Third coordinate array (3D).
  - `nhalo`: Halo width used by node/cell domains.

# Keywords
  - `backend`: Storage backend. Default: `CPU()`.
  - `diff_backend`: Differentiation backend. Default: `AutoForwardDiff()`.
  - `T`: Grid floating-point type.
  - `Tcore`: Deprecated alias for `T` behavior; if provided, overrides numeric type.
  - `compute_metrics`: Enable initial metric computation. Default: `true`.
    Set `compute_metrics=false` with `cache_mode=:off` to disable metric
    allocation/caching entirely.
  - `halo_coords_included`: Whether inputs include halo nodes. Default: `false`.
  - `coordinate_system`: Coordinate-system trait. Default: `CurvilinearCS()`.
  - `basis`: Basis trait. Default: `CartesianBasis()`.
  - `interpolation`: Interpolation mode. Must be `:linear`.
  - `cache_mode`: Metric cache mode (`:eager`, `:lazy`, `:off`). Default: `:eager`.
  - `conserved_metric_scheme`: Conserved face interpolation scheme trait (`EdgeInterpolationOrder1()`, `EdgeInterpolationOrder2()`, `EdgeInterpolationOrder3()`). Default: `EdgeInterpolationOrder3()`.

# Returns
A `DiscreteGrid{N,T,...}` instance with initialized coordinates and metric
cache storage.
"""
function DiscreteGrid(
  x::AbstractVector{TX},
  nhalo::Integer;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T::Type=TX,
  Tcore::Union{Nothing,Type}=nothing,
  compute_metrics::Bool=true,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
) where {TX}
  _validate_discrete_interpolation(interpolation)

  number_type = isnothing(Tcore) ? T : Tcore

  @assert nhalo >= 0

  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  x_nodes = number_type.(x_nodes)

  x_itp = _linear_interpolant(x_nodes)
  x_map = InterpolantMapping{1,typeof(x_itp)}(x_itp)

  return _new_discrete_grid(
    Val(1),
    (; x1=x_map),
    (; x1=x_itp),
    (length(x_nodes) - 1,),
    nhalo;
    backend=backend,
    diff_backend=diff_backend,
    t=zero(number_type),
    params=(;),
    T=number_type,
    compute_metrics=compute_metrics,
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
    conserved_metric_scheme=conserved_metric_scheme,
  )
end

function DiscreteGrid(
  x::AbstractArray{TX,2},
  y::AbstractArray{TX,2},
  nhalo::Integer;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T::Type=TX,
  Tcore::Union{Nothing,Type}=nothing,
  compute_metrics::Bool=true,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
) where {TX}
  size(x) == size(y) ||
    throw(ArgumentError("x and y arrays must have matching dimensions."))
  _validate_discrete_interpolation(interpolation)

  number_type = isnothing(Tcore) ? T : Tcore

  @assert nhalo >= 0

  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  y_nodes = halo_coords_included ? _strip_halo_nodes(y, nhalo) : y

  x_nodes = number_type.(x_nodes)
  y_nodes = number_type.(y_nodes)

  x_itp = _linear_interpolant(x_nodes)
  y_itp = _linear_interpolant(y_nodes)

  x_map = InterpolantMapping{2,typeof(x_itp)}(x_itp)
  y_map = InterpolantMapping{2,typeof(y_itp)}(y_itp)

  return _new_discrete_grid(
    Val(2),
    (; x1=x_map, x2=y_map),
    (; x1=x_itp, x2=y_itp),
    Tuple(size(x_nodes) .- 1),
    nhalo;
    backend=backend,
    diff_backend=diff_backend,
    t=zero(number_type),
    params=(;),
    T=number_type,
    compute_metrics=compute_metrics,
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
    conserved_metric_scheme=conserved_metric_scheme,
  )
end

function DiscreteGrid(
  x::AbstractArray{TX,3},
  y::AbstractArray{TX,3},
  z::AbstractArray{TX,3},
  nhalo::Integer;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T::Type=TX,
  Tcore::Union{Nothing,Type}=nothing,
  compute_metrics::Bool=true,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
  conserved_metric_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
) where {TX}
  (size(x) == size(y) && size(y) == size(z)) ||
    throw(ArgumentError("x, y, and z arrays must have matching dimensions."))
  _validate_discrete_interpolation(interpolation)

  number_type = isnothing(Tcore) ? T : Tcore

  @assert nhalo >= 0

  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  y_nodes = halo_coords_included ? _strip_halo_nodes(y, nhalo) : y
  z_nodes = halo_coords_included ? _strip_halo_nodes(z, nhalo) : z

  x_nodes = number_type.(x_nodes)
  y_nodes = number_type.(y_nodes)
  z_nodes = number_type.(z_nodes)

  x_itp = _linear_interpolant(x_nodes)
  y_itp = _linear_interpolant(y_nodes)
  z_itp = _linear_interpolant(z_nodes)

  x_map = InterpolantMapping{3,typeof(x_itp)}(x_itp)
  y_map = InterpolantMapping{3,typeof(y_itp)}(y_itp)
  z_map = InterpolantMapping{3,typeof(z_itp)}(z_itp)

  return _new_discrete_grid(
    Val(3),
    (; x1=x_map, x2=y_map, x3=z_map),
    (; x1=x_itp, x2=y_itp, x3=z_itp),
    Tuple(size(x_nodes) .- 1),
    nhalo;
    backend=backend,
    diff_backend=diff_backend,
    t=zero(number_type),
    params=(;),
    T=number_type,
    compute_metrics=compute_metrics,
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
    conserved_metric_scheme=conserved_metric_scheme,
  )
end

"""
    update!(grid::DiscreteGrid{N}, t::Real, params::NamedTuple) where {N}

Update discrete grid coordinates and invalidate metric caches for a new state.

# Arguments
  - `grid`: Target discrete grid.
  - `t`: New time value.
  - `params`: New parameter tuple used by mapping functions.

# Returns
`nothing`.
"""
function update!(grid::DiscreteGrid{N}, t::Real, params::NamedTuple) where {N}
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

  grid.state[] = (; t, params)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end

function update!(grid::DiscreteGrid, t::Real=zero(Float64))
  state = grid.state[]
  params = state isa NamedTuple && haskey(state, :params) ? state.params : (;)
  return update!(grid, t, params)
end
