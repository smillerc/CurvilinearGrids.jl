#
# DiscreteGrid
#

using Interpolations

mutable struct DiscreteGrid{C,CS<:CoordinateSystemTrait,BT<:BasisTrait,I} <:
               AbstractMappedOrDiscreteGrid
  core::C
  coordinate_system_trait::CS
  basis_vector_trait::BT
  interpolation::Symbol
  interpolants::I
  state::Any
  metric_caches::UnifiedMetricCaches
end

function _strip_halo_nodes(A, nhalo::Int)
  if nhalo <= 0
    return A
  end
  inds = ntuple(d -> (firstindex(A, d) + nhalo):(lastindex(A, d) - nhalo), ndims(A))
  return @view A[inds...]
end

_linear_interpolant(A) = extrapolate(interpolate(A, BSpline(Linear())), Line())

function _recompute_discrete_cell_metrics!(grid::DiscreteGrid)
  state = grid.state
  t = state.t
  params = state.params

  compute_node_coordinates!(grid.core, t, params)
  compute_centroid_coordinates!(grid.core, t, params)
  compute_cell_metrics!(grid.core, t, params)
  return nothing
end

function _recompute_discrete_face_metrics!(grid::DiscreteGrid)
  state = grid.state
  t = state.t
  params = state.params

  _recompute_discrete_cell_metrics!(grid)
  compute_edge_metrics!(grid.core, t, params)
  return nothing
end

function _refresh_cell_metrics!(grid::DiscreteGrid; include_halo_region::Bool=false)
  _recompute_discrete_cell_metrics!(grid)
  data = getproperty(grid.core, :cell_center_metrics)

  if grid.metric_caches.cell.mode === :off
    return data
  end

  grid.metric_caches.cell.data = deepcopy(data)
  grid.metric_caches.cell.valid = true
  return grid.metric_caches.cell.data
end

function _refresh_face_metrics!(grid::DiscreteGrid; include_halo_region::Bool=false)
  _recompute_discrete_face_metrics!(grid)
  data = getproperty(grid.core, :edge_metrics)

  if grid.metric_caches.face.mode === :off
    return data
  end

  grid.metric_caches.face.data = deepcopy(data)
  grid.metric_caches.face.valid = true
  return grid.metric_caches.face.data
end

function DiscreteGrid(
  core::Union{
    AbstractContinuousCurvilinearGrid1D,
    AbstractContinuousCurvilinearGrid2D,
    AbstractContinuousCurvilinearGrid3D,
  },
  interpolants;
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  interpolation::Symbol=:linear,
  state=(; t=zero(Float64), params=(;)),
  cache_mode::Symbol=:eager,
)
  if interpolation !== :linear
    throw(
      ArgumentError(
        "`DiscreteGrid` currently supports only linear interpolation (`interpolation=:linear`).",
      ),
    )
  end

  caches = _new_metric_caches(cache_mode)
  grid = DiscreteGrid(core, coordinate_system, basis, interpolation, interpolants, state, caches)

  if cache_mode === :eager
    _refresh_cell_metrics!(grid)
    _refresh_face_metrics!(grid)
  end

  return grid
end

function DiscreteGrid(
  x::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  Tcore=T,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {T}
  _, _, _, nhalo, _ = get_gradient_discretization_scheme(discretization_scheme)
  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x

  x_itp = _linear_interpolant(x_nodes)
  x_map(_t, ξ, _p) = x_itp(ξ)
  params = (;)
  celldims = (length(x_nodes) - 1,)

  core = ContinuousCurvilinearGrid1D(
    x_map,
    params,
    celldims,
    discretization_scheme,
    backend,
    diff_backend,
    zero(Float64),
    Tcore;
    compute_metrics=false,
  )

  DiscreteGrid(
    core,
    (; x=x_itp);
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
  )
end

function DiscreteGrid(
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  Tcore=T,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {T}
  size(x) == size(y) || throw(ArgumentError("x and y arrays must have matching dimensions."))

  _, _, _, nhalo, _ = get_gradient_discretization_scheme(discretization_scheme)
  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  y_nodes = halo_coords_included ? _strip_halo_nodes(y, nhalo) : y

  x_itp = _linear_interpolant(x_nodes)
  y_itp = _linear_interpolant(y_nodes)

  x_map(_t, ξ, η, _p) = x_itp(ξ, η)
  y_map(_t, ξ, η, _p) = y_itp(ξ, η)
  params = (;)
  celldims = size(x_nodes) .- 1

  core = ContinuousCurvilinearGrid2D(
    x_map,
    y_map,
    params,
    Tuple(celldims),
    discretization_scheme,
    backend,
    diff_backend,
    zero(Float64),
    Tcore;
    compute_metrics=false,
  )

  DiscreteGrid(
    core,
    (; x=x_itp, y=y_itp);
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
  )
end

function DiscreteGrid(
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  discretization_scheme::Symbol;
  backend=CPU(),
  diff_backend=AutoForwardDiff(),
  Tcore=T,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {T}
  (size(x) == size(y) && size(y) == size(z)) ||
    throw(ArgumentError("x, y, and z arrays must have matching dimensions."))

  _, _, _, nhalo, _ = get_gradient_discretization_scheme(discretization_scheme)
  x_nodes = halo_coords_included ? _strip_halo_nodes(x, nhalo) : x
  y_nodes = halo_coords_included ? _strip_halo_nodes(y, nhalo) : y
  z_nodes = halo_coords_included ? _strip_halo_nodes(z, nhalo) : z

  x_itp = _linear_interpolant(x_nodes)
  y_itp = _linear_interpolant(y_nodes)
  z_itp = _linear_interpolant(z_nodes)

  x_map(_t, ξ, η, ζ, _p) = x_itp(ξ, η, ζ)
  y_map(_t, ξ, η, ζ, _p) = y_itp(ξ, η, ζ)
  z_map(_t, ξ, η, ζ, _p) = z_itp(ξ, η, ζ)
  params = (;)
  celldims = size(x_nodes) .- 1

  core = ContinuousCurvilinearGrid3D(
    x_map,
    y_map,
    z_map,
    params,
    Tuple(celldims),
    discretization_scheme,
    backend,
    diff_backend,
    zero(Float64),
    Tcore;
    compute_metrics=false,
  )

  DiscreteGrid(
    core,
    (; x=x_itp, y=y_itp, z=z_itp);
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
  )
end

function update!(grid::DiscreteGrid, t::Real=zero(Float64))
  params = grid.state.params
  compute_node_coordinates!(grid.core, t, params)
  compute_centroid_coordinates!(grid.core, t, params)
  grid.state = (; t, params)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end

function update!(grid::DiscreteGrid, t::Real, params::NamedTuple)
  compute_node_coordinates!(grid.core, t, params)
  compute_centroid_coordinates!(grid.core, t, params)
  grid.state = (; t, params)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end
