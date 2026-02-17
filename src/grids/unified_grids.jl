"""
Unified grid model (Phase 1, additive):
- MappedGrid
- DiscreteGrid
- OrthogonalGrid
"""

abstract type AbstractUnifiedGrid end
abstract type AbstractMappedOrDiscreteGrid <: AbstractUnifiedGrid end

#
# Traits
#

abstract type CoordinateSystemTrait end
struct CartesianCS <: CoordinateSystemTrait end
struct CylindricalCS <: CoordinateSystemTrait end
struct SphericalCS <: CoordinateSystemTrait end
struct AxisymmetricCS{Axis} <: CoordinateSystemTrait end
struct CurvilinearCS <: CoordinateSystemTrait end

abstract type BasisTrait end
struct CartesianBasis <: BasisTrait end
struct ContravariantBasis <: BasisTrait end
struct CovariantBasis <: BasisTrait end
struct SphericalBasis <: BasisTrait end

#
# Independent metric caches
#

mutable struct UnifiedMetricCache
  data::Any
  valid::Bool
  mode::Symbol
end

mutable struct UnifiedMetricCaches
  cell::UnifiedMetricCache
  face::UnifiedMetricCache
end

function _check_cache_mode(cache_mode::Symbol)
  if cache_mode ∉ (:eager, :lazy, :off)
    throw(
      ArgumentError(
        "Invalid cache mode `$cache_mode`. Expected one of `:eager`, `:lazy`, `:off`.",
      ),
    )
  end
  return cache_mode
end

function _new_metric_caches(cache_mode::Symbol)
  mode = _check_cache_mode(cache_mode)
  UnifiedMetricCaches(
    UnifiedMetricCache(nothing, false, mode), UnifiedMetricCache(nothing, false, mode)
  )
end

#
# Unified grid types
#

mutable struct MappedGrid{L,CS<:CoordinateSystemTrait,BT<:BasisTrait} <:
               AbstractMappedOrDiscreteGrid
  legacy::L
  coordinate_system_trait::CS
  basis_vector_trait::BT
  state::Any
  metric_caches::UnifiedMetricCaches
end

mutable struct DiscreteGrid{L,CS<:CoordinateSystemTrait,BT<:BasisTrait} <:
               AbstractMappedOrDiscreteGrid
  legacy::L
  coordinate_system_trait::CS
  basis_vector_trait::BT
  interpolation::Symbol
  metric_caches::UnifiedMetricCaches
end

mutable struct OrthogonalGrid{L,CS<:CoordinateSystemTrait} <: AbstractUnifiedGrid
  legacy::L
  coordinate_system_trait::CS
  geometry_cache::Any
end

#
# Trait helpers
#

coordinate_system(grid::AbstractUnifiedGrid) = grid.coordinate_system_trait
basis_trait(grid::Union{MappedGrid,DiscreteGrid}) = grid.basis_vector_trait
function basis_trait(::OrthogonalGrid)
  throw(ArgumentError("`basis_trait` is undefined for `OrthogonalGrid`."))
end

coordinate_system(::Type{<:MappedGrid{L,CS}}) where {L,CS} = CS()
coordinate_system(::Type{<:DiscreteGrid{L,CS}}) where {L,CS} = CS()
coordinate_system(::Type{<:OrthogonalGrid{L,CS}}) where {L,CS} = CS()

basis_trait(::Type{<:MappedGrid{L,CS,BT}}) where {L,CS,BT} = BT()
basis_trait(::Type{<:DiscreteGrid{L,CS,BT}}) where {L,CS,BT} = BT()
function basis_trait(::Type{<:OrthogonalGrid})
  throw(ArgumentError("`basis_trait` is undefined for `OrthogonalGrid`."))
end

#
# Legacy trait inference
#

_coordinate_system_from_legacy(::CartesianOrthogonalGrid1D) = CartesianCS()
_coordinate_system_from_legacy(::CylindricalOrthogonalGrid1D) = CylindricalCS()
_coordinate_system_from_legacy(::SphericalOrthogonalGrid1D) = SphericalCS()
_coordinate_system_from_legacy(::SphericalGrid3D) = SphericalCS()
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

_coordinate_system_from_legacy(::AxisymmetricOrthogonalGrid2D) = AxisymmetricCS{:y}()
_coordinate_system_from_legacy(::AbstractCurvilinearGrid) = CurvilinearCS()

_basis_trait_from_legacy(::SphericalBasisCurvilinearGrid3D) = SphericalBasis()
_basis_trait_from_legacy(::AbstractCurvilinearGrid) = ContravariantBasis()

function _is_orthogonal_legacy(mesh)
  mesh isa Union{
    CartesianOrthogonalGrid1D,
    CylindricalOrthogonalGrid1D,
    SphericalOrthogonalGrid1D,
    AxisymmetricOrthogonalGrid2D,
    SphericalGrid3D,
  }
end

#
# Metric cache operations
#

_legacy_cell_metrics(mesh) = getproperty(mesh, :cell_center_metrics)
_legacy_face_metrics(mesh) = getproperty(mesh, :edge_metrics)

function invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)
  grid.metric_caches.cell.valid = false
  grid.metric_caches.cell.data = nothing
  return nothing
end

function invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)
  grid.metric_caches.face.valid = false
  grid.metric_caches.face.data = nothing
  return nothing
end

function refresh_cell_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false)
  if grid.metric_caches.cell.mode === :off
    return _legacy_cell_metrics(grid.legacy)
  end

  grid.metric_caches.cell.data = _legacy_cell_metrics(grid.legacy)
  grid.metric_caches.cell.valid = true
  return grid.metric_caches.cell.data
end

function refresh_face_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false)
  if grid.metric_caches.face.mode === :off
    return _legacy_face_metrics(grid.legacy)
  end

  grid.metric_caches.face.data = _legacy_face_metrics(grid.legacy)
  grid.metric_caches.face.valid = true
  return grid.metric_caches.face.data
end

function cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  if refresh || !grid.metric_caches.cell.valid
    return refresh_cell_metrics!(grid)
  end
  return grid.metric_caches.cell.data
end

function face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  if refresh || !grid.metric_caches.face.valid
    return refresh_face_metrics!(grid)
  end
  return grid.metric_caches.face.data
end

function cell_metrics(::OrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`cell_metrics` is undefined for `OrthogonalGrid`."))
end

function face_metrics(::OrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`face_metrics` is undefined for `OrthogonalGrid`."))
end

#
# Unified constructors from legacy grids
#

function MappedGrid(
  legacy::Union{
    AbstractContinuousCurvilinearGrid1D,
    AbstractContinuousCurvilinearGrid2D,
    AbstractContinuousCurvilinearGrid3D,
  };
  coordinate_system::CoordinateSystemTrait=_coordinate_system_from_legacy(legacy),
  basis::BasisTrait=_basis_trait_from_legacy(legacy),
  state=nothing,
  cache_mode::Symbol=:eager,
)
  caches = _new_metric_caches(cache_mode)
  grid = MappedGrid(legacy, coordinate_system, basis, state, caches)
  if cache_mode === :eager
    refresh_cell_metrics!(grid)
    refresh_face_metrics!(grid)
  end
  return grid
end

function DiscreteGrid(
  legacy::AbstractCurvilinearGrid;
  coordinate_system::CoordinateSystemTrait=_coordinate_system_from_legacy(legacy),
  basis::BasisTrait=_basis_trait_from_legacy(legacy),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
)
  if legacy isa Union{
    AbstractContinuousCurvilinearGrid1D,
    AbstractContinuousCurvilinearGrid2D,
    AbstractContinuousCurvilinearGrid3D,
  }
    throw(ArgumentError("Continuous mapping grids must use `MappedGrid(...)`."))
  end

  if _is_orthogonal_legacy(legacy)
    throw(ArgumentError("Orthogonal legacy grids must use `OrthogonalGrid(...)`."))
  end

  if interpolation !== :linear
    throw(
      ArgumentError(
        "`DiscreteGrid` currently supports only linear interpolation (`interpolation=:linear`).",
      ),
    )
  end

  caches = _new_metric_caches(cache_mode)
  grid = DiscreteGrid(legacy, coordinate_system, basis, interpolation, caches)
  if cache_mode === :eager
    refresh_cell_metrics!(grid)
    refresh_face_metrics!(grid)
  end
  return grid
end

function OrthogonalGrid(
  legacy::Union{
    CartesianOrthogonalGrid1D,
    CylindricalOrthogonalGrid1D,
    SphericalOrthogonalGrid1D,
    AxisymmetricOrthogonalGrid2D,
    SphericalGrid3D,
  };
  coordinate_system::CoordinateSystemTrait=_coordinate_system_from_legacy(legacy),
  geometry_cache=nothing,
)
  OrthogonalGrid(legacy, coordinate_system, geometry_cache)
end

#
# New constructors (phase 1 wrappers around existing constructors)
#

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
  legacy = ContinuousCurvilinearGrid1D(
    x,
    params,
    celldims,
    discretization_scheme,
    backend,
    diff_backend,
    t,
    T;
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
  )
  MappedGrid(
    legacy;
    coordinate_system=coordinate_system,
    basis=basis,
    state=(; t, params),
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
  T=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  legacy = ContinuousCurvilinearGrid2D(
    x,
    y,
    params,
    celldims,
    discretization_scheme,
    backend,
    diff_backend,
    t,
    T;
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
  )
  MappedGrid(
    legacy;
    coordinate_system=coordinate_system,
    basis=basis,
    state=(; t, params),
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
  T=Float64,
  compute_metrics=true,
  global_cell_indices=nothing,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  cache_mode::Symbol=:eager,
)
  legacy = ContinuousCurvilinearGrid3D(
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
    compute_metrics=compute_metrics,
    global_cell_indices=global_cell_indices,
  )
  MappedGrid(
    legacy;
    coordinate_system=coordinate_system,
    basis=basis,
    state=(; t, params),
    cache_mode=cache_mode,
  )
end

function DiscreteGrid(
  x::AbstractVector{T},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=false,
  empty_metrics=false,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {T}
  legacy = CurvilinearGrid1D(
    x,
    discretization_scheme;
    backend=backend,
    is_static=is_static,
    empty_metrics=empty_metrics,
    halo_coords_included=halo_coords_included,
  )
  DiscreteGrid(
    legacy;
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
  is_static=false,
  init_metrics=true,
  empty_metrics=false,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {T}
  legacy = CurvilinearGrid2D(
    x,
    y,
    discretization_scheme;
    backend=backend,
    is_static=is_static,
    init_metrics=init_metrics,
    empty_metrics=empty_metrics,
    halo_coords_included=halo_coords_included,
  )
  DiscreteGrid(
    legacy;
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
  is_static=false,
  init_metrics=true,
  empty_metrics=false,
  halo_coords_included=false,
  coordinate_system::CoordinateSystemTrait=CurvilinearCS(),
  basis::BasisTrait=ContravariantBasis(),
  interpolation::Symbol=:linear,
  cache_mode::Symbol=:eager,
) where {T}
  legacy = CurvilinearGrid3D(
    x,
    y,
    z,
    discretization_scheme;
    backend=backend,
    is_static=is_static,
    init_metrics=init_metrics,
    empty_metrics=empty_metrics,
    halo_coords_included=halo_coords_included,
  )
  DiscreteGrid(
    legacy;
    coordinate_system=coordinate_system,
    basis=basis,
    interpolation=interpolation,
    cache_mode=cache_mode,
  )
end

function OrthogonalGrid(
  x::AbstractVector{T},
  nhalo::Int;
  backend=CPU(),
  coordinate_system::CoordinateSystemTrait=CartesianCS(),
  halo_coords_included=false,
) where {T<:Real}
  legacy = if coordinate_system isa CartesianCS
    CartesianOrthogonalGrid1D(x, nhalo, backend; halo_coords_included=halo_coords_included)
  elseif coordinate_system isa CylindricalCS
    CylindricalOrthogonalGrid1D(x, nhalo, backend; halo_coords_included=halo_coords_included)
  elseif coordinate_system isa SphericalCS
    SphericalOrthogonalGrid1D(x, nhalo, backend; halo_coords_included=halo_coords_included)
  else
    throw(
      ArgumentError(
        "Unsupported coordinate system for 1D orthogonal grid: $(typeof(coordinate_system))",
      ),
    )
  end

  OrthogonalGrid(legacy; coordinate_system=coordinate_system)
end

function OrthogonalGrid(
  r::AbstractVector{T},
  z::AbstractVector{T},
  nhalo::Int;
  backend=CPU(),
  coordinate_system::CoordinateSystemTrait=AxisymmetricCS{:y}(),
  halo_coords_included=false,
) where {T<:Real}
  if !(coordinate_system isa AxisymmetricCS)
    throw(ArgumentError("2D orthogonal constructor currently supports only axisymmetric CS."))
  end
  legacy = AxisymmetricOrthogonalGrid2D(
    r, z, nhalo, backend; halo_coords_included=halo_coords_included
  )
  OrthogonalGrid(legacy; coordinate_system=coordinate_system)
end

function OrthogonalGrid(
  r::AbstractVector{T},
  theta::AbstractVector{T},
  phi::AbstractVector{T},
  nhalo::Int;
  backend=CPU(),
  coordinate_system::CoordinateSystemTrait=SphericalCS(),
  halo_coords_included=false,
) where {T<:Real}
  if !(coordinate_system isa SphericalCS)
    throw(ArgumentError("3D orthogonal constructor currently supports only spherical CS."))
  end
  legacy = SphericalGrid3D(r, theta, phi, nhalo, backend; halo_coords_included=halo_coords_included)
  OrthogonalGrid(legacy; coordinate_system=coordinate_system)
end

#
# Adapter access and delegated geometry API
#

legacy_grid(grid::AbstractUnifiedGrid) = grid.legacy

coords(grid::AbstractUnifiedGrid) = coords(grid.legacy)
coord(grid::AbstractUnifiedGrid, idx) = coord(grid.legacy, idx)
centroids(grid::AbstractUnifiedGrid) = centroids(grid.legacy)
centroid(grid::AbstractUnifiedGrid, idx) = centroid(grid.legacy, idx)
cellvolume(grid::AbstractUnifiedGrid, idx) = cellvolume(grid.legacy, idx)
cellvolumes(grid::AbstractUnifiedGrid) = cellvolumes(grid.legacy)
cellsize(grid::AbstractUnifiedGrid) = cellsize(grid.legacy)
cellsize_withhalo(grid::AbstractUnifiedGrid) = cellsize_withhalo(grid.legacy)

function jacobian_matrix(grid::AbstractMappedOrDiscreteGrid, idx)
  jacobian_matrix(grid.legacy, idx)
end

forward_cell_metrics(grid::AbstractMappedOrDiscreteGrid, idx) = forward_cell_metrics(grid.legacy, idx)
inverse_cell_metrics(grid::AbstractMappedOrDiscreteGrid, idx) = inverse_cell_metrics(grid.legacy, idx)

#
# Updates
#

function _update_continuous_legacy!(
  mesh::Union{
    AbstractContinuousCurvilinearGrid1D,
    AbstractContinuousCurvilinearGrid2D,
    AbstractContinuousCurvilinearGrid3D,
  },
  t,
  params,
)
  compute_node_coordinates!(mesh, t, params)
  compute_centroid_coordinates!(mesh, t, params)
  compute_cell_metrics!(mesh, t, params)
  compute_edge_metrics!(mesh, t, params)
  return nothing
end

function update!(grid::MappedGrid, t, params)
  if grid.legacy isa Union{
    AbstractContinuousCurvilinearGrid1D,
    AbstractContinuousCurvilinearGrid2D,
    AbstractContinuousCurvilinearGrid3D,
  }
    _update_continuous_legacy!(grid.legacy, t, params)
  else
    update!(grid.legacy, t, params)
  end
  grid.state = (; t, params)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end

function update!(grid::MappedGrid, args...; kwargs...)
  update!(grid.legacy, args...; kwargs...)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end

function update!(grid::DiscreteGrid, args...; kwargs...)
  update!(grid.legacy, args...; kwargs...)
  invalidate_cell_metrics!(grid)
  invalidate_face_metrics!(grid)
  return nothing
end
