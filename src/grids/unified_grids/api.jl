#
# Cache API
#

#
# Trait helpers
#

"""
    coordinate_system(grid::AbstractUnifiedGrid)

Return the coordinate-system trait associated with a unified grid.

# Arguments
  - `grid`: Unified grid instance.

# Returns
Coordinate-system trait instance (for example `CurvilinearCS()` or `SphericalCS()`).
"""
coordinate_system(grid::AbstractUnifiedGrid) = coordinate_system(typeof(grid))

"""
    basis_trait(grid::Union{MappedGrid,DiscreteGrid})

Return the basis trait associated with a mapped or discrete unified grid.

# Arguments
  - `grid`: Mapped or discrete unified grid instance.

# Returns
Basis trait instance (`CartesianBasis()` or `SphericalBasis()`).
"""
basis_trait(grid::Union{MappedGrid,DiscreteGrid}) = basis_trait(typeof(grid))
function basis_trait(::OrthogonalGrid)
  throw(ArgumentError("`basis_trait` is undefined for `OrthogonalGrid`."))
end

coordinate_system(::Type{<:MappedGrid{N,T,CS}}) where {N,T,CS} = CS()
coordinate_system(::Type{<:DiscreteGrid{N,T,CS}}) where {N,T,CS} = CS()
coordinate_system(::Type{<:OrthogonalGrid{N,T,L,CS}}) where {N,T,L,CS} = CS()

basis_trait(::Type{<:MappedGrid{N,T,CS,BT}}) where {N,T,CS,BT} = BT()
basis_trait(::Type{<:DiscreteGrid{N,T,CS,BT}}) where {N,T,CS,BT} = BT()
function basis_trait(::Type{<:OrthogonalGrid})
  throw(ArgumentError("`basis_trait` is undefined for `OrthogonalGrid`."))
end

Base.eltype(::MappedGrid{N,T}) where {N,T} = T
Base.eltype(::DiscreteGrid{N,T}) where {N,T} = T
Base.eltype(::OrthogonalGrid{N,T}) where {N,T} = T

"""
    invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)

Mark cell-metric cache entries as stale.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Returns
`nothing`.
"""
function invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)
  if !_has_metric_storage(grid)
    return nothing
  end
  grid.metric_caches.cell.valid = false
  return nothing
end

"""
    invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)

Mark face-metric cache entries as stale.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Returns
`nothing`.
"""
function invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)
  if !_has_metric_storage(grid)
    return nothing
  end
  grid.metric_caches.face.valid = false
  return nothing
end

"""
    refresh_cell_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region=false)

Recompute and return cell metrics for a mapped/discrete unified grid.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `include_halo_region`: Reserved compatibility flag. Default: `false`.

# Returns
Cell metric cache payload.
"""
function refresh_cell_metrics!(
  grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false
)
  _require_metric_storage(grid, "refresh_cell_metrics!")
  _refresh_cell_metrics!(grid; include_halo_region=include_halo_region)
end

"""
    refresh_face_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region=false)

Recompute and return face metrics for a mapped/discrete unified grid.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `include_halo_region`: Reserved compatibility flag. Default: `false`.

# Returns
Face metric cache payload.
"""
function refresh_face_metrics!(
  grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false
)
  _require_metric_storage(grid, "refresh_face_metrics!")
  _refresh_face_metrics!(grid; include_halo_region=include_halo_region)
end

"""
    cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh=false)

Access cell metric cache data.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `refresh`: Recompute before returning data. Default: `false`.

# Returns
Cell metric cache payload.
"""
function cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  _require_metric_storage(grid, "cell_metrics")
  if refresh || !grid.metric_caches.cell.valid
    return refresh_cell_metrics!(grid)
  end
  return grid.metric_caches.cell.data
end

"""
    face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh=false)

Access face metric cache data.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `refresh`: Recompute before returning data. Default: `false`.

# Returns
Face metric cache payload.
"""
function face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool=false)
  _require_metric_storage(grid, "face_metrics")
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
# Geometry API
#

legacy_grid(grid::OrthogonalGrid) = grid.legacy

function coords(grid::Union{MappedGrid{1},DiscreteGrid{1}})
  @views grid.node_coordinates[1][grid.iterators.node.domain]
end

function coords(grid::Union{MappedGrid{2},DiscreteGrid{2}})
  @views (
    grid.node_coordinates[1][grid.iterators.node.domain],
    grid.node_coordinates[2][grid.iterators.node.domain],
  )
end

function coords(grid::Union{MappedGrid{3},DiscreteGrid{3}})
  @views (
    grid.node_coordinates[1][grid.iterators.node.domain],
    grid.node_coordinates[2][grid.iterators.node.domain],
    grid.node_coordinates[3][grid.iterators.node.domain],
  )
end

coords(grid::OrthogonalGrid) = coords(grid.legacy)

function coord(grid::Union{MappedGrid{1},DiscreteGrid{1}}, (i,)::NTuple{1,Int})
  @SVector [grid.node_coordinates[1][i]]
end
function coord(grid::Union{MappedGrid{2},DiscreteGrid{2}}, (i, j)::NTuple{2,Int})
  @SVector [grid.node_coordinates[1][i, j], grid.node_coordinates[2][i, j]]
end
function coord(grid::Union{MappedGrid{3},DiscreteGrid{3}}, (i, j, k)::NTuple{3,Int})
  @SVector [
    grid.node_coordinates[1][i, j, k],
    grid.node_coordinates[2][i, j, k],
    grid.node_coordinates[3][i, j, k],
  ]
end

coord(grid::OrthogonalGrid, idx) = coord(grid.legacy, idx)

function centroids(grid::Union{MappedGrid{1},DiscreteGrid{1}})
  @views grid.centroid_coordinates[1][grid.iterators.cell.domain]
end

function centroids(grid::Union{MappedGrid{2},DiscreteGrid{2}})
  @views (
    grid.centroid_coordinates[1][grid.iterators.cell.domain],
    grid.centroid_coordinates[2][grid.iterators.cell.domain],
  )
end

function centroids(grid::Union{MappedGrid{3},DiscreteGrid{3}})
  @views (
    grid.centroid_coordinates[1][grid.iterators.cell.domain],
    grid.centroid_coordinates[2][grid.iterators.cell.domain],
    grid.centroid_coordinates[3][grid.iterators.cell.domain],
  )
end

centroids(grid::OrthogonalGrid) = centroids(grid.legacy)

function centroid(grid::Union{MappedGrid{1},DiscreteGrid{1}}, (i,)::NTuple{1,Int})
  @SVector [grid.centroid_coordinates[1][i]]
end
function centroid(grid::Union{MappedGrid{2},DiscreteGrid{2}}, (i, j)::NTuple{2,Int})
  @SVector [grid.centroid_coordinates[1][i, j], grid.centroid_coordinates[2][i, j]]
end
function centroid(grid::Union{MappedGrid{3},DiscreteGrid{3}}, (i, j, k)::NTuple{3,Int})
  @SVector [
    grid.centroid_coordinates[1][i, j, k],
    grid.centroid_coordinates[2][i, j, k],
    grid.centroid_coordinates[3][i, j, k],
  ]
end

centroid(grid::OrthogonalGrid, idx) = centroid(grid.legacy, idx)

@inline _promote_real_tuple(idx::Tuple{Vararg{Real,N}}) where {N} = promote(idx...)

@inline function _state_for_eval(grid::AbstractMappedOrDiscreteGrid)
  state = grid.state[]
  if !(state isa NamedTuple && haskey(state, :t) && haskey(state, :params))
    throw(ArgumentError("Grid state is missing `(t, params)` and cannot be evaluated."))
  end
  return state.t, state.params
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::Tuple{Vararg{Real,1}}
)
  ξ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [grid.mapping_functions.x1(t, ξ..., params)]
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::Tuple{Vararg{Real,2}}
)
  ξη = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [
    grid.mapping_functions.x1(t, ξη..., params), grid.mapping_functions.x2(t, ξη..., params)
  ]
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{3},DiscreteGrid{3}}, idx::Tuple{Vararg{Real,3}}
)
  ξηζ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [
    grid.mapping_functions.x1(t, ξηζ..., params),
    grid.mapping_functions.x2(t, ξηζ..., params),
    grid.mapping_functions.x3(t, ξηζ..., params),
  ]
end

@inline function _continuous_forward_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _require_metric_functions(grid, "forward_cell_metrics")
  ξηζ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  F = _as_smatrix(Val(N), grid.metric_functions_cache.forward.jacobian(t, ξηζ..., params))
  return SMatrix{N,N,T,N * N}(Tuple(F))
end

@inline function _continuous_inverse_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  Jinv = inv(_continuous_forward_jacobian(grid, idx))
  return SMatrix{N,N,T,N * N}(Tuple(Jinv))
end

function coord(
  grid::Union{MappedGrid{N},DiscreteGrid{N}}, idx::Tuple{Vararg{Real,N}}
) where {N}
  _continuous_coord(grid, idx)
end

@inline function _cell_forward_metric_at(grid::AbstractMappedOrDiscreteGrid, idx)
  cm = cell_metrics(grid)
  cm.forward[idx...]
end

@inline function _cell_inverse_metric_at(grid::AbstractMappedOrDiscreteGrid, idx)
  cm = cell_metrics(grid)
  cm.inverse[idx...]
end

@inline function _cell_forward_metric_at(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  F = _continuous_forward_jacobian(grid, idx)
  return Metric(F, det(F))
end

@inline function _cell_inverse_metric_at(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  G = _continuous_inverse_jacobian(grid, idx)
  return Metric(G, det(G))
end

"""
    forward_cell_metrics(grid, idx)

Return the forward cell metric payload at `idx` as a `Metric`.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: `CartesianIndex`, integer tuple (discrete coordinate), or real tuple
    (continuous coordinate).
"""
@inline function forward_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex
)
  forward_cell_metrics(grid, idx.I)
end
@inline function forward_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::NTuple{N,Int}
) where {N}
  _cell_forward_metric_at(grid, idx)
end
@inline function forward_cell_metrics(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _cell_forward_metric_at(grid, idx)
end

"""
    inverse_cell_metrics(grid, idx)

Return the inverse cell metric payload at `idx` as a `Metric`.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: `CartesianIndex`, integer tuple (discrete coordinate), or real tuple
    (continuous coordinate).
"""
@inline function inverse_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex
)
  inverse_cell_metrics(grid, idx.I)
end
@inline function inverse_cell_metrics(
  grid::Union{MappedGrid,DiscreteGrid}, idx::NTuple{N,Int}
) where {N}
  _cell_inverse_metric_at(grid, idx)
end
@inline function inverse_cell_metrics(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _cell_inverse_metric_at(grid, idx)
end

function forward_cell_metrics(::OrthogonalGrid, idx)
  throw(ArgumentError("`forward_cell_metrics` is undefined for `OrthogonalGrid`."))
end

function inverse_cell_metrics(::OrthogonalGrid, idx)
  throw(ArgumentError("`inverse_cell_metrics` is undefined for `OrthogonalGrid`."))
end

"""
    cellvolume(grid, idx)

Compute cell volume at a given index.

# Arguments
  - `grid`: Unified grid instance.
  - `idx`: `CartesianIndex` or tuple index.
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

cellvolume(grid::OrthogonalGrid, idx) = cellvolume(grid.legacy, idx)

@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  return _cell_forward_metric_at(grid, idx).J
end
@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  return T(det(_continuous_forward_jacobian(grid, idx)))
end

@inline _radial_centroid_1d(grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::NTuple{1,Int}) = grid.centroid_coordinates[1][idx...]
@inline _radial_centroid_1d(grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::Tuple{Vararg{Real,1}}) = _continuous_coord(
  grid, idx
)[1]

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
  r = _radial_centroid_1d(grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  r = _radial_centroid_1d(grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::NTuple{2,Int},
)
  r = centroid(grid, idx)[1]
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  r = _continuous_coord(grid, idx)[1]
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::NTuple{2,Int},
) where {Axis}
  r = _axisymmetric_radius(cs, grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
) where {Axis}
  r = _axisymmetric_radius(cs, grid, idx)
  return (2π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{1,Int}
)
  r = _radial_centroid_1d(grid, idx)
  return (4π * r^2) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  r = _radial_centroid_1d(grid, idx)
  return (4π * r^2) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::NTuple{3,Int},
)
  c = centroid(grid, idx)
  r = c[1]
  θ = c[2]
  return (r^2 * sin(θ)) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::Tuple{Vararg{Real,3}},
)
  c = _continuous_coord(grid, idx)
  r = c[1]
  θ = c[2]
  return (r^2 * sin(θ)) * _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, ::AbstractMappedOrDiscreteGrid, idx::NTuple{2,Int}
)
  throw(ArgumentError("Spherical cellvolume dispatch is not implemented for 2D grids."))
end

function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  ::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
)
  throw(ArgumentError("Spherical cellvolume dispatch is not implemented for 2D grids."))
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

function cellvolumes(grid::Union{MappedGrid,DiscreteGrid})
  volumes = zeros(eltype(grid), size(grid.iterators.cell.domain))

  for (idx0, idx1) in zip(CartesianIndices(volumes), grid.iterators.cell.domain)
    volumes[idx0] = cellvolume(grid, idx1)
  end

  return volumes
end

cellvolumes(grid::OrthogonalGrid) = cellvolumes(grid.legacy)

cellsize(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.domain)
cellsize(grid::OrthogonalGrid) = cellsize(grid.legacy)

cellsize_withhalo(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.full)
cellsize_withhalo(grid::OrthogonalGrid) = cellsize_withhalo(grid.legacy)

"""
    jacobian_matrix(grid, idx)

Return the forward Jacobian matrix at a given cell index.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: `CartesianIndex`, integer tuple, or real tuple.
"""
function jacobian_matrix(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex)
  jacobian_matrix(grid, idx.I)
end
function jacobian_matrix(grid::Union{MappedGrid,DiscreteGrid}, idx::NTuple{N,Int}) where {N}
  _cell_forward_metric_at(grid, idx).jacobian_matrix
end
function jacobian_matrix(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _continuous_forward_jacobian(grid, idx)
end
