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
  - `discretization_scheme`: Discretization scheme object.
  - `discretization_scheme_name`: Symbol name of the discretization scheme.
  - `iterators`: Node/cell iterator bundle for local/global indexing.
  - `interpolation`: Interpolation mode symbol (`:linear`).
  - `interpolants`: Cached interpolation objects.
  - `state`: Mutable reference to runtime state (`t`, `params`).
  - `metric_caches`: Independent cell and face metric caches.
"""
struct DiscreteGrid{
  N,
  T,
  CS<:CoordinateSystemTrait,
  BT<:BasisTrait,
  IP,
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
  )

  return nothing
end

function _recompute_discrete_face_metrics!(
  grid::DiscreteGrid{N,T}; include_halo_region::Bool=false
) where {N,T}
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
  )

  return nothing
end

function _refresh_cell_metrics!(grid::DiscreteGrid; include_halo_region::Bool=false)
  _recompute_discrete_cell_metrics!(grid; include_halo_region=include_halo_region)
  data = grid.metric_caches.cell.data

  if grid.metric_caches.cell.mode === :off
    return data
  end

  grid.metric_caches.cell.valid = true
  return data
end

function _refresh_face_metrics!(grid::DiscreteGrid; include_halo_region::Bool=false)
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
  discretization_scheme::Symbol;
  backend,
  diff_backend,
  t,
  params,
  T::Type,
  coordinate_system::CoordinateSystemTrait,
  basis::BasisTrait,
  interpolation::Symbol,
  cache_mode::Symbol,
) where {N}
  _check_unified_basis_trait(basis)
  _validate_discrete_interpolation(interpolation)

  components = _build_unified_components(
    Val(N), mapping_functions, celldims, discretization_scheme, backend, diff_backend, T
  )

  caches = _new_metric_caches(
    cache_mode, components.cell_metric_storage, components.face_metric_storage
  )
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

  if cache_mode === :eager
    _refresh_cell_metrics!(grid)
    _refresh_face_metrics!(grid)
  end

  return grid
end

"""
    DiscreteGrid(x[, y[, z]], discretization_scheme; kwargs...)

Construct a discrete unified grid from coordinate arrays using linear
interpolation.

# Arguments
  - `x`: First coordinate array.
  - `y`: Second coordinate array (2D/3D).
  - `z`: Third coordinate array (3D).
  - `discretization_scheme`: Gradient scheme symbol (for example `:meg6`).

# Keywords
  - `backend`: Storage backend. Default: `CPU()`.
  - `diff_backend`: Differentiation backend. Default: `AutoForwardDiff()`.
  - `T`: Grid floating-point type.
  - `Tcore`: Deprecated alias for `T` behavior; if provided, overrides numeric type.
  - `halo_coords_included`: Whether inputs include halo nodes. Default: `false`.
  - `coordinate_system`: Coordinate-system trait. Default: `CurvilinearCS()`.
  - `basis`: Basis trait. Default: `CartesianBasis()`.
  - `interpolation`: Interpolation mode. Must be `:linear`.
  - `cache_mode`: Metric cache mode (`:eager`, `:lazy`, `:off`). Default: `:eager`.

# Returns
A `DiscreteGrid{N,T,...}` instance with initialized coordinates and metric
cache storage.
"""
function DiscreteGrid(
  x::AbstractVector{TX},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T::Type=TX,
  Tcore::Union{Nothing,Type}=nothing,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {TX}
  _validate_discrete_interpolation(interpolation)

  number_type = isnothing(Tcore) ? T : Tcore

  _, _, _, nhalo, _ = get_gradient_discretization_scheme(discretization_scheme)
  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  x_nodes = number_type.(x_nodes)

  x_itp = _linear_interpolant(x_nodes)
  x_map(_t, ξ, _p) = x_itp(ξ)

  return _new_discrete_grid(
    Val(1),
    (; x=x_map),
    (; x=x_itp),
    (length(x_nodes) - 1,),
    discretization_scheme;
    backend=backend,
    diff_backend=diff_backend,
    t=zero(number_type),
    params=(;),
    T=number_type,
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
  )
end

function DiscreteGrid(
  x::AbstractArray{TX,2},
  y::AbstractArray{TX,2},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T::Type=TX,
  Tcore::Union{Nothing,Type}=nothing,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {TX}
  size(x) == size(y) ||
    throw(ArgumentError("x and y arrays must have matching dimensions."))
  _validate_discrete_interpolation(interpolation)

  number_type = isnothing(Tcore) ? T : Tcore

  _, _, _, nhalo, _ = get_gradient_discretization_scheme(discretization_scheme)
  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  y_nodes = halo_coords_included ? _strip_halo_nodes(y, nhalo) : y

  x_nodes = number_type.(x_nodes)
  y_nodes = number_type.(y_nodes)

  x_itp = _linear_interpolant(x_nodes)
  y_itp = _linear_interpolant(y_nodes)

  x_map(_t, ξ, η, _p) = x_itp(ξ, η)
  y_map(_t, ξ, η, _p) = y_itp(ξ, η)

  return _new_discrete_grid(
    Val(2),
    (; x=x_map, y=y_map),
    (; x=x_itp, y=y_itp),
    Tuple(size(x_nodes) .- 1),
    discretization_scheme;
    backend=backend,
    diff_backend=diff_backend,
    t=zero(number_type),
    params=(;),
    T=number_type,
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
  )
end

function DiscreteGrid(
  x::AbstractArray{TX,3},
  y::AbstractArray{TX,3},
  z::AbstractArray{TX,3},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  T::Type=TX,
  Tcore::Union{Nothing,Type}=nothing,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=CartesianBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {TX}
  (size(x) == size(y) && size(y) == size(z)) ||
    throw(ArgumentError("x, y, and z arrays must have matching dimensions."))
  _validate_discrete_interpolation(interpolation)

  number_type = isnothing(Tcore) ? T : Tcore

  _, _, _, nhalo, _ = get_gradient_discretization_scheme(discretization_scheme)
  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  y_nodes = halo_coords_included ? _strip_halo_nodes(y, nhalo) : y
  z_nodes = halo_coords_included ? _strip_halo_nodes(z, nhalo) : z

  x_nodes = number_type.(x_nodes)
  y_nodes = number_type.(y_nodes)
  z_nodes = number_type.(z_nodes)

  x_itp = _linear_interpolant(x_nodes)
  y_itp = _linear_interpolant(y_nodes)
  z_itp = _linear_interpolant(z_nodes)

  x_map(_t, ξ, η, ζ, _p) = x_itp(ξ, η, ζ)
  y_map(_t, ξ, η, ζ, _p) = y_itp(ξ, η, ζ)
  z_map(_t, ξ, η, ζ, _p) = z_itp(ξ, η, ζ)

  return _new_discrete_grid(
    Val(3),
    (; x=x_map, y=y_map, z=z_map),
    (; x=x_itp, y=y_itp, z=z_itp),
    Tuple(size(x_nodes) .- 1),
    discretization_scheme;
    backend=backend,
    diff_backend=diff_backend,
    t=zero(number_type),
    params=(;),
    T=number_type,
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
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

function update!(grid::DiscreteGrid, t::Real=zero(Float64))
  state = grid.state[]
  params = state isa NamedTuple && haskey(state, :params) ? state.params : (;)
  return update!(grid, t, params)
end
