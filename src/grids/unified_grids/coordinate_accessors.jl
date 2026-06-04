#
# Coordinate accessors
#

#
# Coordinate arrays and plotting coordinates
#

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

"""
    cartesian_coordinates(grid::AbstractOrthogonalGrid)

Return orthogonal-grid node coordinates in Cartesian form.

For coordinate systems already stored in Cartesian-like form, this returns
`coords(grid)` unchanged. For spherical 3D grids, this converts `(r, θ, ϕ)` to
`(x, y, z)`.
"""
@inline function cartesian_coordinates(grid::AbstractOrthogonalGrid)
  _cartesian_coordinates(coordinate_system(grid), coords(grid))
end

"""
    cartesian_centroids(grid)

Return 2-D cell-center coordinates in Cartesian plotting form.

For Cartesian-like orthogonal grids, tensor-product centroid vectors are expanded
to matrices. For spherical 2-D grids, `(r, θ)` centroids are converted to
Cartesian meridional coordinates `(x, z)`. For mapped/discrete grids, cached
centroid arrays are returned in Cartesian form when needed.
"""
@inline function cartesian_centroids(
  grid::Union{MappedGrid{2},DiscreteGrid{2},OrthogonalGrid{2}}
)
  return _cartesian_centroids(coordinate_system(grid), centroids(grid))
end

@inline _cartesian_coordinates(::CoordinateSystemTrait, q) = q
@inline _cartesian_coordinates(::CartesianCS, q) = q
@inline _cartesian_coordinates(::CurvilinearCS, q) = q
@inline _cartesian_coordinates(::AxisymmetricCS, q) = q
@inline _cartesian_coordinates(::CylindricalCS, q) = q

@inline function _tensor_product_coordinates(
  x::AbstractVector{Tx}, y::AbstractVector{Ty}
) where {Tx,Ty}
  xx = [x[i] for i in eachindex(x), _ in eachindex(y)]
  yy = [y[j] for _ in eachindex(x), j in eachindex(y)]
  return xx, yy
end

@inline _tensor_product_coordinates(x::AbstractArray, y::AbstractArray) = (x, y)

@inline _cartesian_centroids(::CartesianCS, q::NTuple{2,Any}) =
  _tensor_product_coordinates(q...)
@inline _cartesian_centroids(::CurvilinearCS, q::NTuple{2,Any}) =
  _tensor_product_coordinates(q...)
@inline _cartesian_centroids(::AxisymmetricCS, q::NTuple{2,Any}) =
  _tensor_product_coordinates(q...)
@inline _cartesian_centroids(::CylindricalCS, q::NTuple{2,Any}) =
  _tensor_product_coordinates(q...)
@inline _cartesian_centroids(::CoordinateSystemTrait, q::NTuple{2,Any}) = q

@inline function _cartesian_coordinates(
  ::SphericalCS, q::NTuple{2,<:AbstractVector{T}}
) where {T}
  r, θ = q
  R = reshape(r, :, 1)
  sinθ = reshape(sin.(θ), 1, :)
  cosθ = reshape(cos.(θ), 1, :)
  x = R .* sinθ
  z = R .* cosθ
  return x, z
end

@inline function _cartesian_centroids(
  ::SphericalCS, q::NTuple{2,<:AbstractVector{T}}
) where {T}
  r, θ = q
  R = reshape(r, :, 1)
  sinθ = reshape(sin.(θ), 1, :)
  cosθ = reshape(cos.(θ), 1, :)
  return R .* sinθ, R .* cosθ
end

@inline function _cartesian_centroids(::SphericalCS, q::NTuple{2,<:AbstractArray})
  r, θ = q
  return @. r * sin(θ), @. r * cos(θ)
end

@inline function _cartesian_coordinates(
  ::SphericalCS, q::NTuple{3,<:AbstractVector{T}}
) where {T}
  r, θ, ϕ = q
  R = reshape(r, :, 1, 1)
  sinθ = reshape(sin.(θ), 1, :, 1)
  cosθ = reshape(cos.(θ), 1, :, 1)
  sinϕ = reshape(sin.(ϕ), 1, 1, :)
  cosϕ = reshape(cos.(ϕ), 1, 1, :)
  oneϕ = reshape(one.(ϕ), 1, 1, :)

  x = R .* sinθ .* cosϕ
  y = R .* sinθ .* sinϕ
  z = R .* cosθ .* oneϕ

  return x, y, z
end

"""
    coord(grid, idx)

Return the native physical coordinate at a node or continuous computational
coordinate.

For `MappedGrid` and `DiscreteGrid`, integer tuple and `CartesianIndex` inputs
address stored node coordinates. Real-valued tuples evaluate the underlying
mapping at that computational coordinate, for example `coord(grid, (5.25, 7.5))`.
Orthogonal-grid methods return coordinates assembled from the native node
coordinate arrays.

See also [`coords`](@ref), [`centroid`](@ref), and [`jacobian_matrix`](@ref).
"""

#
# Discrete coordinate accessors
#

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

"""
    centroids(grid)

Return views of the cell-center coordinates over the non-halo interior domain.

The returned coordinates are expressed in the grid's native physical coordinate
system:
- Cartesian-like grids return Cartesian components,
- orthogonal spherical grids return `(r, θ)` or `(r, θ, ϕ)`,
- axisymmetric grids return their native meridional coordinates.

For `MappedGrid` and `DiscreteGrid`, the values come from the centroid cache.
For `OrthogonalGrid`, the corresponding methods are defined in
`GridTypes.jl` and follow the same native-coordinate convention.
"""
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

"""
    face_coordinates(grid)

Return the stored high-side face-center coordinate arrays for mapped/discrete
grids.

The layout is `face_coordinates(grid)[axis][component][I]`. The index `I`
matches `face_metrics(grid)[axis].conserved[I]`; low-side faces are read from
the neighboring high-side face index.
"""
@inline face_coordinates(grid::Union{MappedGrid,DiscreteGrid}) = grid.face_coordinates

"""
    centroid(grid, idx)

Return the native physical coordinate of the cell center at `idx`.

For mapped/discrete grids, `idx` is an integer cell index. Orthogonal-grid
methods are also available and return `SVector`s assembled from the native
centroid coordinate arrays.

See also [`centroids`](@ref), [`face_coordinate`](@ref), and
[`basis_transfer_matrix`](@ref).
"""
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
#
# Continuous coordinate evaluation
#

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
