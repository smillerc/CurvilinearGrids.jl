#
# Cache API
#

#
# Trait helpers
#

coordinate_system(grid::AbstractUnifiedGrid) = coordinate_system(typeof(grid))
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

function invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)
  grid.metric_caches.cell.valid = false
  return nothing
end

function invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)
  grid.metric_caches.face.valid = false
  return nothing
end

function refresh_cell_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false)
  _refresh_cell_metrics!(grid; include_halo_region=include_halo_region)
end

function refresh_face_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool=false)
  _refresh_face_metrics!(grid; include_halo_region=include_halo_region)
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
# Geometry API
#

legacy_grid(grid::OrthogonalGrid) = grid.legacy

coords(grid::Union{MappedGrid{1},DiscreteGrid{1}}) =
  @views grid.node_coordinates[1][grid.iterators.node.domain]

coords(grid::Union{MappedGrid{2},DiscreteGrid{2}}) = @views (
  grid.node_coordinates[1][grid.iterators.node.domain],
  grid.node_coordinates[2][grid.iterators.node.domain],
)

coords(grid::Union{MappedGrid{3},DiscreteGrid{3}}) = @views (
  grid.node_coordinates[1][grid.iterators.node.domain],
  grid.node_coordinates[2][grid.iterators.node.domain],
  grid.node_coordinates[3][grid.iterators.node.domain],
)

coords(grid::OrthogonalGrid) = coords(grid.legacy)

coord(grid::Union{MappedGrid{1},DiscreteGrid{1}}, (i,)::NTuple{1,Int}) =
  @SVector [grid.node_coordinates[1][i]]
coord(grid::Union{MappedGrid{2},DiscreteGrid{2}}, (i, j)::NTuple{2,Int}) =
  @SVector [grid.node_coordinates[1][i, j], grid.node_coordinates[2][i, j]]
coord(grid::Union{MappedGrid{3},DiscreteGrid{3}}, (i, j, k)::NTuple{3,Int}) = @SVector [
  grid.node_coordinates[1][i, j, k],
  grid.node_coordinates[2][i, j, k],
  grid.node_coordinates[3][i, j, k],
]

coord(grid::OrthogonalGrid, idx) = coord(grid.legacy, idx)

centroids(grid::Union{MappedGrid{1},DiscreteGrid{1}}) =
  @views grid.centroid_coordinates[1][grid.iterators.cell.domain]

centroids(grid::Union{MappedGrid{2},DiscreteGrid{2}}) = @views (
  grid.centroid_coordinates[1][grid.iterators.cell.domain],
  grid.centroid_coordinates[2][grid.iterators.cell.domain],
)

centroids(grid::Union{MappedGrid{3},DiscreteGrid{3}}) = @views (
  grid.centroid_coordinates[1][grid.iterators.cell.domain],
  grid.centroid_coordinates[2][grid.iterators.cell.domain],
  grid.centroid_coordinates[3][grid.iterators.cell.domain],
)

centroids(grid::OrthogonalGrid) = centroids(grid.legacy)

centroid(grid::Union{MappedGrid{1},DiscreteGrid{1}}, (i,)::NTuple{1,Int}) =
  @SVector [grid.centroid_coordinates[1][i]]
centroid(grid::Union{MappedGrid{2},DiscreteGrid{2}}, (i, j)::NTuple{2,Int}) =
  @SVector [grid.centroid_coordinates[1][i, j], grid.centroid_coordinates[2][i, j]]
centroid(grid::Union{MappedGrid{3},DiscreteGrid{3}}, (i, j, k)::NTuple{3,Int}) = @SVector [
  grid.centroid_coordinates[1][i, j, k],
  grid.centroid_coordinates[2][i, j, k],
  grid.centroid_coordinates[3][i, j, k],
]

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
  return @SVector [grid.mapping_functions.x(t, ξ..., params)]
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::Tuple{Vararg{Real,2}}
)
  ξη = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [
    grid.mapping_functions.x(t, ξη..., params),
    grid.mapping_functions.y(t, ξη..., params),
  ]
end

@inline function _continuous_coord(
  grid::Union{MappedGrid{3},DiscreteGrid{3}}, idx::Tuple{Vararg{Real,3}}
)
  ξηζ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  return @SVector [
    grid.mapping_functions.x(t, ξηζ..., params),
    grid.mapping_functions.y(t, ξηζ..., params),
    grid.mapping_functions.z(t, ξηζ..., params),
  ]
end

@inline function _continuous_forward_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  ξηζ = _promote_real_tuple(idx)
  t, params = _state_for_eval(grid)
  F = _as_smatrix(Val(N), grid.metric_functions_cache.forward.jacobian(t, ξηζ..., params))
  return SMatrix{N,N,T,N * N}(Tuple(F))
end

@inline function _continuous_inverse_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  F = _continuous_forward_jacobian(grid, idx)
  G = inv(F)
  return SMatrix{N,N,T,N * N}(Tuple(G))
end

coord(
  grid::Union{MappedGrid{N},DiscreteGrid{N}}, idx::Tuple{Vararg{Real,N}}
) where {N} = _continuous_coord(grid, idx)

centroid(
  grid::Union{MappedGrid{N},DiscreteGrid{N}}, idx::Tuple{Vararg{Real,N}}
) where {N} = _continuous_coord(grid, idx)

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
  return _metric_from_jacobian(F, T(det(F)))
end

@inline function _cell_inverse_metric_at(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  F = _continuous_forward_jacobian(grid, idx)
  G = _continuous_inverse_jacobian(grid, idx)
  return _metric_from_jacobian(G, T(det(F)))
end

cellvolume(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex) = cellvolume(grid, idx.I)
cellvolume(
  grid::Union{MappedGrid{N,T,CS,BT},DiscreteGrid{N,T,CS,BT}}, idx::NTuple{N,Int}
) where {N,T,CS,BT} = _cellvolume_dispatch(CS(), BT(), grid, idx)
cellvolume(
  grid::Union{MappedGrid{N,T,CS,BT},DiscreteGrid{N,T,CS,BT}}, idx::Tuple{Vararg{Real,N}}
) where {N,T,CS,BT} = _cellvolume_dispatch(CS(), BT(), grid, _promote_real_tuple(idx))

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

@inline _radial_centroid_1d(grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::NTuple{1,Int}) =
  grid.centroid_coordinates[1][idx...]
@inline _radial_centroid_1d(
  grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::Tuple{Vararg{Real,1}}
) = _continuous_coord(grid, idx)[1]

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:x}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::NTuple{2,Int}
)
  return grid.centroid_coordinates[2][idx...]
end
@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:x}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::Tuple{Vararg{Real,2}}
)
  return _continuous_coord(grid, idx)[2]
end

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:y}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::NTuple{2,Int}
)
  return grid.centroid_coordinates[1][idx...]
end
@inline function _axisymmetric_radius(
  ::AxisymmetricCS{:y}, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::Tuple{Vararg{Real,2}}
)
  return _continuous_coord(grid, idx)[1]
end

@inline function _axisymmetric_radius(
  ::AxisymmetricCS{Axis}, grid::AbstractMappedOrDiscreteGrid, idx::Tuple{Vararg{Real,2}}
) where {Axis}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis`. Supported values are `:x` and `:y`.",
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
  return (2 * π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS, ::CartesianBasis, grid::AbstractMappedOrDiscreteGrid, idx::Tuple{Vararg{Real,1}}
)
  r = _radial_centroid_1d(grid, idx)
  return (2 * π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS, ::CartesianBasis, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::NTuple{2,Int}
)
  r = centroid(grid, idx)[1]
  return (2 * π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS, ::CartesianBasis, grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx::Tuple{Vararg{Real,2}}
)
  r = _continuous_coord(grid, idx)[1]
  return (2 * π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::NTuple{2,Int},
) where {Axis}
  r = _axisymmetric_radius(cs, grid, idx)
  return (2 * π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
) where {Axis}
  r = _axisymmetric_radius(cs, grid, idx)
  return (2 * π * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{1,Int}
)
  r = _radial_centroid_1d(grid, idx)
  return (4 * π * r * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::Tuple{Vararg{Real,1}}
)
  r = _radial_centroid_1d(grid, idx)
  return (4 * π * r * r) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::Union{MappedGrid{3},DiscreteGrid{3}}, idx::NTuple{3,Int}
)
  c = centroid(grid, idx)
  r = c[1]
  θ = c[2]
  return (r * r * sin(θ)) * _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::Union{MappedGrid{3},DiscreteGrid{3}}, idx::Tuple{Vararg{Real,3}}
)
  c = _continuous_coord(grid, idx)
  r = c[1]
  θ = c[2]
  return (r * r * sin(θ)) * _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, ::AbstractMappedOrDiscreteGrid, idx::NTuple{2,Int}
)
  throw(
    ArgumentError(
      "Spherical cellvolume dispatch is not implemented for 2D grids.",
    ),
  )
end

function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, ::AbstractMappedOrDiscreteGrid, idx::Tuple{Vararg{Real,2}}
)
  throw(
    ArgumentError(
      "Spherical cellvolume dispatch is not implemented for 2D grids.",
    ),
  )
end

function _cellvolume_dispatch(
  cs::CoordinateSystemTrait,
  bt::SphericalBasis,
  ::AbstractMappedOrDiscreteGrid,
  ::NTuple{N,Int},
) where {N}
  throw(
    ArgumentError(
      "Unsupported basis/coordinate-system combination: $(typeof(bt)) with $(typeof(cs)).",
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
      "Unsupported basis/coordinate-system combination: $(typeof(bt)) with $(typeof(cs)).",
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

jacobian_matrix(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex) = jacobian_matrix(
  grid, idx.I
)
function jacobian_matrix(grid::Union{MappedGrid,DiscreteGrid}, idx::NTuple{N,Int}) where {N}
  _cell_forward_metric_at(grid, idx).jacobian_matrix
end
function jacobian_matrix(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  _continuous_forward_jacobian(grid, idx)
end

function forward_cell_metrics(grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx)
  F = _cell_forward_metric_at(grid, idx isa CartesianIndex ? idx.I : idx).jacobian_matrix
  return (; x=(; ξ=F[1, 1]),)
end

function forward_cell_metrics(grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx)
  F = _cell_forward_metric_at(grid, idx isa CartesianIndex ? idx.I : idx).jacobian_matrix
  return (; x=(; ξ=F[1, 1], η=F[1, 2]), y=(; ξ=F[2, 1], η=F[2, 2]))
end

function forward_cell_metrics(grid::Union{MappedGrid{3},DiscreteGrid{3}}, idx)
  F = _cell_forward_metric_at(grid, idx isa CartesianIndex ? idx.I : idx).jacobian_matrix
  return (
    ;
    x=(; ξ=F[1, 1], η=F[1, 2], ζ=F[1, 3]),
    y=(; ξ=F[2, 1], η=F[2, 2], ζ=F[2, 3]),
    z=(; ξ=F[3, 1], η=F[3, 2], ζ=F[3, 3]),
  )
end

function inverse_cell_metrics(grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx)
  m = _cell_inverse_metric_at(grid, idx isa CartesianIndex ? idx.I : idx)
  G = m.jacobian_matrix
  Ghat = G .* m.J
  return (; ξ=(; x₁=G[1, 1]), ξ̂=(; x₁=Ghat[1, 1]))
end

function inverse_cell_metrics(grid::Union{MappedGrid{2},DiscreteGrid{2}}, idx)
  m = _cell_inverse_metric_at(grid, idx isa CartesianIndex ? idx.I : idx)
  G = m.jacobian_matrix
  Ghat = G .* m.J
  return (
    ;
    ξ=(; x₁=G[1, 1], x₂=G[1, 2]),
    η=(; x₁=G[2, 1], x₂=G[2, 2]),
    ξ̂=(; x₁=Ghat[1, 1], x₂=Ghat[1, 2]),
    η̂=(; x₁=Ghat[2, 1], x₂=Ghat[2, 2]),
  )
end

function inverse_cell_metrics(grid::Union{MappedGrid{3},DiscreteGrid{3}}, idx)
  m = _cell_inverse_metric_at(grid, idx isa CartesianIndex ? idx.I : idx)
  G = m.jacobian_matrix
  Ghat = G .* m.J
  return (
    ;
    ξ=(; x₁=G[1, 1], x₂=G[1, 2], x₃=G[1, 3]),
    η=(; x₁=G[2, 1], x₂=G[2, 2], x₃=G[2, 3]),
    ζ=(; x₁=G[3, 1], x₂=G[3, 2], x₃=G[3, 3]),
    ξ̂=(; x₁=Ghat[1, 1], x₂=Ghat[1, 2], x₃=Ghat[1, 3]),
    η̂=(; x₁=Ghat[2, 1], x₂=Ghat[2, 2], x₃=Ghat[2, 3]),
    ζ̂=(; x₁=Ghat[3, 1], x₂=Ghat[3, 2], x₃=Ghat[3, 3]),
  )
end
