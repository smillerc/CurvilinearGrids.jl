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
coordinate_system(grid::AbstractOrthogonalGrid) = coordinate_system(typeof(grid))

"""
    basis_trait(grid::Union{MappedGrid,DiscreteGrid})

Return the basis trait associated with a mapped or discrete unified grid.

# Arguments
  - `grid`: Mapped or discrete unified grid instance.

# Returns
Basis trait instance (`CartesianBasis()` or `SphericalBasis()`).
"""
basis_trait(grid::Union{MappedGrid,DiscreteGrid}) = basis_trait(typeof(grid))
function basis_trait(::AbstractOrthogonalGrid)
  throw(ArgumentError("`basis_trait` is undefined for `AbstractOrthogonalGrid`."))
end

coordinate_system(::Type{<:MappedGrid{N,T,CS}}) where {N,T,CS} = CS()
coordinate_system(::Type{<:DiscreteGrid{N,T,CS}}) where {N,T,CS} = CS()
coordinate_system(::Type{<:OrthogonalGrid{N,T,CS}}) where {N,T,CS} = CS()

basis_trait(::Type{<:MappedGrid{N,T,CS,BT}}) where {N,T,CS,BT} = BT()
basis_trait(::Type{<:DiscreteGrid{N,T,CS,BT}}) where {N,T,CS,BT} = BT()
function basis_trait(::Type{<:AbstractOrthogonalGrid})
  throw(ArgumentError("`basis_trait` is undefined for `AbstractOrthogonalGrid`."))
end

"""
    basis_transfer_matrix(grid, from_q, to_q)

Return the local basis-transfer matrix that maps vector components stored at
`from_q` into the local basis at `to_q`.

For mapped and discrete grids, the transfer uses the grid's
`coordinate_system(grid)` and `basis_trait(grid)`. For orthogonal grids, the
transfer is implied by the coordinate system because `basis_trait` is undefined.

Supported public cases in the current API are:
- identity transfer for Cartesian-basis mapped/discrete grids,
- spherical-basis transfer for 2-D and 3-D spherical mapped/discrete grids,
- identity transfer for orthogonal Cartesian, 1-D cylindrical, 1-D spherical,
  and 2-D axisymmetric meridional grids,
- spherical transfer for 2-D and 3-D orthogonal spherical grids.

Unsupported basis/coordinate-system combinations throw `ArgumentError` rather
than silently falling back.

# Arguments
  - `grid`: Unified grid instance.
  - `from_q`: Native physical coordinates at the donor location.
  - `to_q`: Native physical coordinates at the receiver location.

# Returns
Static square matrix `R_(to<-from)` such that `v_to = R_(to<-from) * v_from`.

See also [`face_coordinate`](@ref), [`centroid`](@ref), and the four-argument
`basis_transfer_matrix(from_grid, from_q, to_grid, to_q)` overload for
inter-grid transfers.
"""
function basis_transfer_matrix end

@inline function basis_transfer_matrix(
  from_grid::Union{AbstractUnifiedGrid,AbstractOrthogonalGrid},
  from_q::Tuple{Vararg{Real,N}},
  to_grid::Union{AbstractUnifiedGrid,AbstractOrthogonalGrid},
  to_q::Tuple{Vararg{Real,N}},
) where {N}
  basis_transfer_matrix(from_grid, SVector{N}(from_q), to_grid, SVector{N}(to_q))
end

@inline function basis_transfer_matrix(
  from_grid::Union{AbstractUnifiedGrid,AbstractOrthogonalGrid},
  from_q::SVector{N,T1},
  to_grid::Union{AbstractUnifiedGrid,AbstractOrthogonalGrid},
  to_q::SVector{N,T2},
) where {N,T1,T2}
  T = promote_type(T1, T2)
  from_qT = SVector{N,T}(from_q)
  to_qT = SVector{N,T}(to_q)
  Qfrom = _grid_basis_to_cartesian_matrix(from_grid, from_qT, Val(N), T)
  Qto = _grid_basis_to_cartesian_matrix(to_grid, to_qT, Val(N), T)
  return SMatrix{N,N,T,N * N}(transpose(Qto) * Qfrom)
end

@inline function basis_transfer_matrix(
  grid::Union{MappedGrid{N},DiscreteGrid{N}},
  from_q::Tuple{Vararg{Real,N}},
  to_q::Tuple{Vararg{Real,N}},
) where {N}
  basis_transfer_matrix(grid, SVector{N}(from_q), SVector{N}(to_q))
end

@inline function basis_transfer_matrix(
  grid::Union{MappedGrid{N},DiscreteGrid{N}}, from_q::SVector{N,T1}, to_q::SVector{N,T2}
) where {N,T1,T2}
  return basis_transfer_matrix(grid, from_q, grid, to_q)
end

@inline function basis_transfer_matrix(
  grid::OrthogonalGrid{N}, from_q::Tuple{Vararg{Real,N}}, to_q::Tuple{Vararg{Real,N}}
) where {N}
  basis_transfer_matrix(grid, SVector{N}(from_q), SVector{N}(to_q))
end

@inline function basis_transfer_matrix(
  grid::OrthogonalGrid{N}, from_q::SVector{N,T1}, to_q::SVector{N,T2}
) where {N,T1,T2}
  return basis_transfer_matrix(grid, from_q, grid, to_q)
end

@inline _identity_basis_transfer(::Val{N}, ::Type{T}) where {N,T} = one(SMatrix{N,N,T})

@inline function _grid_basis_to_cartesian_matrix(
  grid::Union{MappedGrid{N},DiscreteGrid{N}}, q::SVector{N,T}, ::Val{N}, ::Type{T}
) where {N,T}
  return _mapped_basis_to_cartesian_matrix(
    coordinate_system(grid), basis_trait(grid), q, Val(N), T
  )
end

@inline function _grid_basis_to_cartesian_matrix(
  grid::OrthogonalGrid{N}, q::SVector{N,T}, ::Val{N}, ::Type{T}
) where {N,T}
  return _orthogonal_basis_to_cartesian_matrix(coordinate_system(grid), q, Val(N), T)
end

@inline function _mapped_basis_to_cartesian_matrix(
  cs::CoordinateSystemTrait, ::CartesianBasis, q::SVector{N,T}, ::Val{N}, ::Type{T}
) where {N,T}
  _ = cs
  _ = q
  return _identity_basis_transfer(Val(N), T)
end

@inline function _basis_to_cartesian_matrix(
  ::SphericalCS, ::Val{2}, q::SVector{2,T}, ::Type{T}
) where {T}
  _, θ = q
  sθ = sin(θ)
  cθ = cos(θ)
  return @SMatrix [sθ cθ; cθ -sθ]
end

@inline function _basis_to_cartesian_matrix(
  ::SphericalCS, ::Val{3}, q::SVector{3,T}, ::Type{T}
) where {T}
  _, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return @SMatrix [
    sθ*cϕ cθ*cϕ -sϕ
    sθ*sϕ cθ*sϕ cϕ
    cθ -sθ 0
  ]
end

@inline function _mapped_basis_to_cartesian_matrix(
  ::SphericalCS, ::SphericalBasis, q::SVector{2,T}, ::Val{2}, ::Type{T}
) where {T}
  return _basis_to_cartesian_matrix(SphericalCS(), Val(2), q, T)
end

@inline function _mapped_basis_to_cartesian_matrix(
  ::SphericalCS, ::SphericalBasis, q::SVector{3,T}, ::Val{3}, ::Type{T}
) where {T}
  return _basis_to_cartesian_matrix(SphericalCS(), Val(3), q, T)
end

@inline function _mapped_basis_to_cartesian_matrix(
  cs::CoordinateSystemTrait, bt::SphericalBasis, q::SVector{N,T}, ::Val{N}, ::Type{T}
) where {N,T}
  _ = q
  throw(
    ArgumentError(
      "Unsupported public basis transfer for $(typeof(bt)) with coordinate system $(typeof(cs)) and N=$N.",
    ),
  )
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  ::CartesianCS, q::SVector{N,T}, ::Val{N}, ::Type{T}
) where {N,T}
  _ = q
  return _identity_basis_transfer(Val(N), T)
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  ::CurvilinearCS, q::SVector{N,T}, ::Val{N}, ::Type{T}
) where {N,T}
  _ = q
  return _identity_basis_transfer(Val(N), T)
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  ::AxisymmetricCS, q::SVector{2,T}, ::Val{2}, ::Type{T}
) where {T}
  _ = q
  return _identity_basis_transfer(Val(2), T)
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  ::CylindricalCS, q::SVector{1,T}, ::Val{1}, ::Type{T}
) where {T}
  _ = q
  return _identity_basis_transfer(Val(1), T)
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  ::SphericalCS, q::SVector{1,T}, ::Val{1}, ::Type{T}
) where {T}
  _ = q
  return _identity_basis_transfer(Val(1), T)
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  ::SphericalCS, q::SVector{2,T}, ::Val{2}, ::Type{T}
) where {T}
  return _basis_to_cartesian_matrix(SphericalCS(), Val(2), q, T)
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  ::SphericalCS, q::SVector{3,T}, ::Val{3}, ::Type{T}
) where {T}
  return _basis_to_cartesian_matrix(SphericalCS(), Val(3), q, T)
end

@inline function _orthogonal_basis_to_cartesian_matrix(
  cs::CoordinateSystemTrait, q::SVector{N,T}, ::Val{N}, ::Type{T}
) where {N,T}
  _ = q
  throw(
    ArgumentError(
      "Unsupported public basis transfer for orthogonal grid coordinate system $(typeof(cs)) and N=$N.",
    ),
  )
end

Base.eltype(::MappedGrid{N,T}) where {N,T} = T
Base.eltype(::DiscreteGrid{N,T}) where {N,T} = T

# Backward-compatible identity constructor for orthogonal grids.
OrthogonalGrid(grid::OrthogonalGrid) = grid

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

function cell_metrics(::AbstractOrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`cell_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

function face_metrics(::AbstractOrthogonalGrid; refresh::Bool=false)
  throw(ArgumentError("`face_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

#
# Geometry API
#

"""
    FaceFluxGeometry{N,T}

Solver-facing face geometry for mapped or discrete grids.

# Fields
  - `coordinate`: Native physical coordinate of the face center.
  - `metric_vector`: Outward-oriented active conserved face metric vector in the
    grid's physical basis. This is the solver-facing face metric `S^alpha`,
    obtained from the active row of
    `face_metrics(grid)[axis].conserved[idx].jacobian_matrix`.
  - `area`: Magnitude of `metric_vector`.
  - `normal`: `metric_vector / area` when `area > 0`, otherwise the zero vector.
"""
struct FaceFluxGeometry{N,T}
  coordinate::SVector{N,T}
  metric_vector::SVector{N,T}
  area::T
  normal::SVector{N,T}
end

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

"""
    InverseCoordinateResult{N,T}

Result payload returned by inverse coordinate mapping routines.

# Fields
  - `coordinate`: Computational coordinate estimate.
  - `converged`: Whether the nonlinear solve converged.
  - `iterations`: Number of Newton iterations performed.
  - `residual_norm`: Infinity norm of the final physical-space residual.
"""
struct InverseCoordinateResult{N,T}
  coordinate::SVector{N,T}
  converged::Bool
  iterations::Int
  residual_norm::T
end

@inline function _computational_bounds(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}
) where {N,T}
  first_node = first(grid.iterators.global_domain.node.domain).I
  last_node = last(grid.iterators.global_domain.node.domain).I

  ξmin = ntuple(i -> T(first_node[i] - grid.nhalo), N)
  ξmax = ntuple(i -> T(last_node[i] - grid.nhalo), N)
  return SVector{N,T}(ξmin), SVector{N,T}(ξmax)
end

@inline function _physical_bounds(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}
) where {N,T}
  mins = MVector{N,T}(undef)
  maxs = MVector{N,T}(undef)
  @inbounds for n in 1:N
    mins[n] = typemax(T)
    maxs[n] = typemin(T)
  end

  @inbounds for I in grid.iterators.node.domain
    for n in 1:N
      x = grid.node_coordinates[n][I]
      mins[n] = min(mins[n], x)
      maxs[n] = max(maxs[n], x)
    end
  end

  return SVector{N,T}(mins), SVector{N,T}(maxs)
end

@inline function _seed_from_bounding_box(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  x_target::SVector{N,T},
  bounds_tolerance::T,
) where {N,T}
  xmin, xmax = _physical_bounds(grid)
  ξmin, ξmax = _computational_bounds(grid)

  @inbounds for n in 1:N
    lo = xmin[n] - bounds_tolerance
    hi = xmax[n] + bounds_tolerance
    if !(lo <= x_target[n] <= hi)
      throw(
        ArgumentError(
          "Physical point lies outside the grid axis-aligned bounding box: x[$n]=$(x_target[n]) ∉ [$lo, $hi].",
        ),
      )
    end
  end

  ξ0 = ntuple(n -> begin
    span = xmax[n] - xmin[n]
    α = span > eps(T) ? (x_target[n] - xmin[n]) / span : T(0.5)
    ξmin[n] + clamp(α, zero(T), one(T)) * (ξmax[n] - ξmin[n])
  end, N)
  return SVector{N,T}(ξ0)
end

@inline function _normalize_inverse_point(::Type{T}, x::Tuple{Vararg{Real,N}}) where {N,T}
  x_promoted = _promote_real_tuple(x)
  SVector{N,T}(ntuple(i -> T(x_promoted[i]), N))
end

@inline function _normalize_inverse_guess(
  ::Type{T}, guess::Tuple{Vararg{Real,N}}
) where {N,T}
  guess_promoted = _promote_real_tuple(guess)
  SVector{N,T}(ntuple(i -> T(guess_promoted[i]), N))
end

@inline function _inverse_coordinate_result(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  x_phys::Tuple{Vararg{Real,N}};
  guess::Union{Nothing,Tuple{Vararg{Real,N}}}=nothing,
  tol::Real=sqrt(eps(T)),
  maxiters::Int=25,
  damping::Bool=true,
  bounds_tolerance::Real=sqrt(eps(T)),
) where {N,T}
  if maxiters <= 0
    throw(ArgumentError("`maxiters` must be positive."))
  end
  tolT = T(tol)
  if tolT <= zero(T)
    throw(ArgumentError("`tol` must be positive."))
  end
  bounds_tolT = T(bounds_tolerance)
  if bounds_tolT < zero(T)
    throw(ArgumentError("`bounds_tolerance` must be non-negative."))
  end

  x_target = _normalize_inverse_point(T, x_phys)
  ξ = if isnothing(guess)
    _seed_from_bounding_box(grid, x_target, bounds_tolT)
  else
    _normalize_inverse_guess(T, guess)
  end

  residual_norm = typemax(T)
  converged = false
  iters = 0
  for iter in 1:maxiters
    iters = iter
    r = _continuous_coord(grid, Tuple(ξ)) - x_target
    residual_norm = norm(r, Inf)
    if residual_norm <= tolT
      converged = true
      break
    end

    Jxξ = _continuous_forward_jacobian(grid, Tuple(ξ))
    Δξ = Jxξ \ r
    ξtrial = ξ - Δξ

    if damping
      λ = one(T)
      improved = false
      while λ >= T(1) / T(64)
        ξcandidate = ξ - λ * Δξ
        rcandidate = _continuous_coord(grid, Tuple(ξcandidate)) - x_target
        if norm(rcandidate, Inf) < residual_norm
          ξtrial = ξcandidate
          improved = true
          break
        end
        λ *= T(0.5)
      end
      if !improved
        ξtrial = ξ - Δξ
      end
    end

    ξ = ξtrial
  end

  if !converged
    r = _continuous_coord(grid, Tuple(ξ)) - x_target
    residual_norm = norm(r, Inf)
    converged = residual_norm <= tolT
  end

  return InverseCoordinateResult{N,T}(ξ, converged, iters, residual_norm)
end

"""
    computational_coordinate(grid, x_phys; kwargs...)

Solve the inverse map `x(ξ) = x_phys` for mapped or discrete unified grids.

The routine uses Newton iterations on the forward mapping with Jacobian
`∂x/∂ξ`, seeded from an axis-aligned physical-space bounding box.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `x_phys`: Physical-space point as an `NTuple{N,Real}`.

# Keywords
  - `guess`: Optional initial guess for `ξ`. Default: `nothing` (bounding-box seed).
  - `tol`: Infinity-norm residual tolerance. Default: `sqrt(eps(T))`.
  - `maxiters`: Maximum Newton iterations. Default: `25`.
  - `damping`: Enable simple backtracking damping. Default: `true`.
  - `bounds_tolerance`: Tolerance for bounding-box containment check. Default: `sqrt(eps(T))`.
  - `throw_on_failure`: Throw if nonlinear solve does not converge. Default: `true`.
  - `return_result`: Return `InverseCoordinateResult` instead of only `ξ`. Default: `false`.

# Returns
Either computational coordinate `SVector{N,T}` or `InverseCoordinateResult{N,T}`
when `return_result=true`.
"""
@inline function computational_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  x_phys::Tuple{Vararg{Real,N}};
  guess::Union{Nothing,Tuple{Vararg{Real,N}}}=nothing,
  tol::Real=sqrt(eps(T)),
  maxiters::Int=25,
  damping::Bool=true,
  bounds_tolerance::Real=sqrt(eps(T)),
  throw_on_failure::Bool=true,
  return_result::Bool=false,
) where {N,T}
  result = _inverse_coordinate_result(
    grid,
    x_phys;
    guess=guess,
    tol=tol,
    maxiters=maxiters,
    damping=damping,
    bounds_tolerance=bounds_tolerance,
  )
  if throw_on_failure && !result.converged
    throw(
      ErrorException(
        "Inverse mapping failed to converge in $(result.iterations) iterations (residual=$(result.residual_norm)).",
      ),
    )
  end
  return return_result ? result : result.coordinate
end

function computational_coordinate(grid::AbstractOrthogonalGrid, x_phys; kwargs...)
  throw(
    ArgumentError("`computational_coordinate` is undefined for `AbstractOrthogonalGrid`.")
  )
end

@inline function _cell_forward_metric_at(grid::AbstractMappedOrDiscreteGrid, idx)
  cm = cell_metrics(grid)
  cm.forward[idx...]
end

@inline function _cell_inverse_metric_at(grid::AbstractMappedOrDiscreteGrid, idx)
  cm = cell_metrics(grid)
  cm.inverse[idx...]
end

@inline function _cell_center_computational_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  Iglobal = grid.iterators.global_domain.cell.full[CartesianIndex(idx)]
  half = T(0.5)
  return ntuple(d -> T(Iglobal.I[d] - grid.nhalo) + half, N)
end

@inline function _cell_forward_metric_at(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  if !_has_metric_storage(grid)
    return _cell_forward_metric_at(grid, _cell_center_computational_coordinate(grid, idx))
  end
  cm = cell_metrics(grid)
  return cm.forward[idx...]
end

@inline function _cell_inverse_metric_at(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  if !_has_metric_storage(grid)
    return _cell_inverse_metric_at(grid, _cell_center_computational_coordinate(grid, idx))
  end
  cm = cell_metrics(grid)
  return cm.inverse[idx...]
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
  return Metric(G, inv(det(G)))
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

function forward_cell_metrics(::AbstractOrthogonalGrid, idx)
  throw(ArgumentError("`forward_cell_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

function inverse_cell_metrics(::AbstractOrthogonalGrid, idx)
  throw(ArgumentError("`inverse_cell_metrics` is undefined for `AbstractOrthogonalGrid`."))
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

"""
    cell_jacobian(grid, idx)

Return the geometric Jacobian factor used by diffusion operators at cell index
`idx`.

For mapped/discrete grids this is `abs(forward_cell_metrics(grid, idx).J)`.
For orthogonal grids it is the reduced volumetric measure associated with the
coordinate family, for example `r`, `r^2`, or `r^2 sin(θ)`.
"""
@inline function cell_jacobian(grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex)
  cell_jacobian(grid, idx.I)
end
@inline function cell_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  abs(forward_cell_metrics(grid, idx).J)
end

@inline function _orth_cell_coord(
  grid::OrthogonalGrid{N}, dim::Int, idx::NTuple{N,Int}
) where {N}
  grid.centroid_coordinates[dim][idx[dim]]
end
@inline function _orth_face_coord(
  grid::OrthogonalGrid{N}, dim::Int, idx::NTuple{N,Int}
) where {N}
  grid.node_coordinates[dim][idx[dim]]
end

@inline function _axisymmetric_radial_dim(::AxisymmetricCS{:y}, ::Val{2})
  1
end
@inline function _axisymmetric_radial_dim(::AxisymmetricCS{:x}, ::Val{2})
  2
end
function _axisymmetric_radial_dim(::AxisymmetricCS{Axis}, ::Val{2}) where {Axis}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis`. Supported values are `:x` and `:y`."
    ),
  )
end

@inline function cell_jacobian(grid::AbstractOrthogonalGrid, idx::CartesianIndex)
  cell_jacobian(grid, idx.I)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{N,T,CartesianCS}, idx::NTuple{N,Int}
) where {N,T}
  one(T)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{1,T,CylindricalCS}, idx::NTuple{1,Int}
) where {T}
  _orth_cell_coord(grid, 1, idx)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{2,T,AxisymmetricCS{Axis}}, idx::NTuple{2,Int}
) where {T,Axis}
  ridx = _axisymmetric_radial_dim(AxisymmetricCS{Axis}(), Val(2))
  _orth_cell_coord(grid, ridx, idx)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{1,T,SphericalCS}, idx::NTuple{1,Int}
) where {T}
  r = _orth_cell_coord(grid, 1, idx)
  r^2
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{2,T,SphericalCS}, idx::NTuple{2,Int}
) where {T}
  r = _orth_cell_coord(grid, 1, idx)
  θ = _orth_cell_coord(grid, 2, idx)
  r^2 * sin(θ)
end
@inline function cell_jacobian(
  grid::OrthogonalGrid{3,T,SphericalCS}, idx::NTuple{3,Int}
) where {T}
  r = _orth_cell_coord(grid, 1, idx)
  θ = _orth_cell_coord(grid, 2, idx)
  r^2 * sin(θ)
end
function cell_jacobian(grid::AbstractOrthogonalGrid, idx::NTuple{N,Int}) where {N}
  throw(
    ArgumentError(
      "`cell_jacobian` is undefined for $(typeof(grid)) with $(N)-D integer indices."
    ),
  )
end

"""
    face_metric_coefficient(grid, dim, idx)

Return the axis-aligned face coefficient `J g^{dd}` used by diffusion operators
for face axis `dim` at face index `idx`.

For orthogonal grids this is the scalar coefficient that belongs with the face
measure in the native coordinate system. For mapped/discrete grids this is
computed from the cached forward and inverse face metrics.
"""
@inline function face_metric_coefficient(
  grid::Union{MappedGrid,DiscreteGrid}, dim::Int, idx::CartesianIndex
)
  face_metric_coefficient(grid, dim, idx.I)
end
@inline function face_metric_coefficient(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, dim::Int, idx::NTuple{N,Int}
) where {N,T}
  1 <= dim <= N ||
    throw(ArgumentError("face axis dim=$dim is invalid for $N-D unified grid"))
  fm = face_metrics(grid)[dim]
  J = abs(fm.forward[idx...].J)
  G = fm.inverse[idx...].jacobian_matrix
  gdd = zero(T)
  @inbounds for m in 1:N
    gdd += G[dim, m]^2
  end
  return J * gdd
end

@inline function face_metric_coefficient(
  grid::AbstractOrthogonalGrid, dim::Int, idx::CartesianIndex
)
  face_metric_coefficient(grid, dim, idx.I)
end
@inline function face_metric_coefficient(
  grid::OrthogonalGrid{N,T,CartesianCS}, dim::Int, idx::NTuple{N,Int}
) where {N,T}
  1 <= dim <= N ||
    throw(ArgumentError("face axis dim=$dim is invalid for $N-D orthogonal grid"))
  one(T)
end
@inline function face_metric_coefficient(
  grid::OrthogonalGrid{1,T,CylindricalCS}, dim::Int, idx::NTuple{1,Int}
) where {T}
  dim == 1 || throw(ArgumentError("face axis dim=$dim is invalid for 1-D cylindrical grid"))
  _orth_face_coord(grid, 1, idx)
end
@inline function face_metric_coefficient(
  grid::OrthogonalGrid{2,T,AxisymmetricCS{Axis}}, dim::Int, idx::NTuple{2,Int}
) where {T,Axis}
  1 <= dim <= 2 ||
    throw(ArgumentError("face axis dim=$dim is invalid for 2-D axisymmetric grid"))
  ridx = _axisymmetric_radial_dim(AxisymmetricCS{Axis}(), Val(2))
  if dim == ridx
    return _orth_face_coord(grid, ridx, idx)
  end
  return _orth_cell_coord(grid, ridx, idx)
end
@inline function face_metric_coefficient(
  grid::OrthogonalGrid{1,T,SphericalCS}, dim::Int, idx::NTuple{1,Int}
) where {T}
  dim == 1 || throw(ArgumentError("face axis dim=$dim is invalid for 1-D spherical grid"))
  r = _orth_face_coord(grid, 1, idx)
  r^2
end
@inline function face_metric_coefficient(
  grid::OrthogonalGrid{2,T,SphericalCS}, dim::Int, idx::NTuple{2,Int}
) where {T}
  1 <= dim <= 2 ||
    throw(ArgumentError("face axis dim=$dim is invalid for 2-D spherical grid"))
  if dim == 1
    r = _orth_face_coord(grid, 1, idx)
    θ = _orth_cell_coord(grid, 2, idx)
    return r^2 * sin(θ)
  end
  θ = _orth_face_coord(grid, 2, idx)
  return sin(θ)
end
@inline function face_metric_coefficient(
  grid::OrthogonalGrid{3,T,SphericalCS}, dim::Int, idx::NTuple{3,Int}
) where {T}
  1 <= dim <= 3 ||
    throw(ArgumentError("face axis dim=$dim is invalid for 3-D spherical grid"))
  if dim == 1
    r = _orth_face_coord(grid, 1, idx)
    θ = _orth_cell_coord(grid, 2, idx)
    return r^2 * sin(θ)
  elseif dim == 2
    θ = _orth_face_coord(grid, 2, idx)
    return sin(θ)
  end
  θ = _orth_cell_coord(grid, 2, idx)
  s = sin(θ)
  abs(s) > sqrt(eps(T)) ||
    throw(DomainError(θ, "Spherical phi coefficient undefined at sin(theta)=0"))
  return inv(s)
end
function face_metric_coefficient(
  grid::AbstractOrthogonalGrid, dim::Int, idx::NTuple{N,Int}
) where {N}
  throw(
    ArgumentError(
      "`face_metric_coefficient` is undefined for $(typeof(grid)) with $(N)-D integer indices.",
    ),
  )
end

@inline function _face_loc_axis_side(::Val{N}, loc::Symbol) where {N}
  # Fast-path canonical symbols used in hot loops to avoid String/lowercase allocations.
  axis, side = if loc === :ilo
    (1, :lo)
  elseif loc === :ihi
    (1, :hi)
  elseif loc === :jlo
    (2, :lo)
  elseif loc === :jhi
    (2, :hi)
  elseif loc === :klo
    (3, :lo)
  elseif loc === :khi
    (3, :hi)
  else
    throw(
      ArgumentError(
        "Unsupported face location `$loc`. Expected one of `:ilo/:ihi`, `:jlo/:jhi`, or `:klo/:khi`.",
      ),
    )
  end

  if axis > N
    throw(
      ArgumentError(
        "Face location `$loc` maps to axis $axis, which is invalid for $N-D grid."
      ),
    )
  end
  return axis, side
end

@inline function _face_location_symbol(axis::Int, side::Symbol, ::Val{N}) where {N}
  if axis == 1
    return side === :lo ? :ilo : :ihi
  elseif axis == 2 && N >= 2
    return side === :lo ? :jlo : :jhi
  elseif axis == 3 && N >= 3
    return side === :lo ? :klo : :khi
  end
  throw(ArgumentError("Invalid axis/side combination `(axis=$axis, side=$side)` for N=$N."))
end

"""
    face_coordinate(grid, idx, loc)

Return the native physical coordinate at the face center selected by `loc` for
cell index `idx`.

For orthogonal grids, the coordinate on the active face axis comes from the
node location while tangential coordinates remain cell-centered. For
mapped/discrete grids, the coordinate is evaluated from the underlying mapping
at the face midpoint in computational space.

# Arguments
  - `grid`: Unified grid instance.
  - `idx`: Cell index as `CartesianIndex` or `NTuple{N,Int}`.
  - `loc`: Face selector symbol (for example `:ilo`, `:ihi`, `:jlo`, `:jhi`,
    `:klo`, `:khi`).

# Returns
An `SVector` giving the face-center coordinate in the grid's physical
coordinate system.
"""
@inline function face_coordinate(
  grid::AbstractUnifiedGrid, idx::CartesianIndex, loc::Symbol
)
  face_coordinate(grid, idx.I, loc)
end

@inline function face_coordinate(
  grid::OrthogonalGrid{N,T}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  axis, side = _face_loc_axis_side(Val(N), loc)
  face_idx = ntuple(d -> d == axis ? idx[d] + (side === :hi ? 1 : 0) : idx[d], N)
  return SVector{N,T}(
    ntuple(
      d -> d == axis ? _orth_face_coord(grid, d, face_idx) : _orth_cell_coord(grid, d, idx),
      N,
    ),
  )
end

@inline function face_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  axis, side = _face_loc_axis_side(Val(N), loc)
  face_idx = _face_flux_cache_index(idx, axis, side)
  coords = face_coordinates(grid)
  return SVector{N,T}(ntuple(d -> coords[axis][d][face_idx...], Val(N)))
end

@inline function _face_evaluation_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  Icell::CartesianIndex{N},
  axis::Int,
  side::Symbol,
) where {N,T}
  Iglobal = grid.iterators.global_domain.cell.full[Icell]
  half = T(0.5)
  base = ntuple(d -> T(Iglobal.I[d] - grid.nhalo) + half, N)
  offset = side === :hi ? half : -half
  return ntuple(d -> d == axis ? base[d] + offset : base[d], N)
end

@inline function _face_mapped_coordinate(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  Icell::CartesianIndex{N},
  axis::Int,
  side::Symbol,
) where {N,T}
  ξηζ_face = _face_evaluation_coordinate(grid, Icell, axis, side)
  q = _continuous_coord(grid, ξηζ_face)
  return SVector{N,T}(ntuple(d -> T(q[d]), N))
end

@inline function _face_forward_jacobian(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}},
  Icell::CartesianIndex{N},
  axis::Int,
  side::Symbol,
) where {N,T}
  _require_metric_functions(grid, "face_area")
  ξηζ_face = _face_evaluation_coordinate(grid, Icell, axis, side)
  t, params = _state_for_eval(grid)
  Fraw = grid.metric_functions_cache.forward.jacobian(t, ξηζ_face..., params)
  return SMatrix{N,N,T,N * N}(Tuple(Fraw))
end

@inline function _face_coord_to_cartesian(
  ::CoordinateSystemTrait, q::SVector{N,T}
) where {N,T}
  return q
end

@inline function _face_coord_to_cartesian(::CartesianCS, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _face_coord_to_cartesian(::CurvilinearCS, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _face_coord_to_cartesian(::AxisymmetricCS, q::SVector{2,T}) where {T}
  return q
end

@inline function _face_coord_to_cartesian(::CylindricalCS, q::SVector{2,T}) where {T}
  return q
end

@inline function _face_coord_to_cartesian(::SphericalCS, q::SVector{2,T}) where {T}
  r, θ = q
  return SVector{2,T}(r * sin(θ), r * cos(θ))
end

@inline function _face_coord_to_cartesian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return SVector{3,T}(r * sθ * cϕ, r * sθ * sϕ, r * cθ)
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CoordinateSystemTrait, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CartesianCS, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CurvilinearCS, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::AxisymmetricCS, q::SVector{2,T}
) where {T}
  return one(SMatrix{2,2,T})
end

@inline function _face_coord_to_cartesian_jacobian(
  ::CylindricalCS, q::SVector{2,T}
) where {T}
  return one(SMatrix{2,2,T})
end

@inline function _face_coord_to_cartesian_jacobian(::SphericalCS, q::SVector{2,T}) where {T}
  r, θ = q
  sθ = sin(θ)
  cθ = cos(θ)
  return SMatrix{2,2,T,4}(sθ, r * cθ, cθ, -r * sθ)
end

@inline function _face_coord_to_cartesian_jacobian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return SMatrix{3,3,T,9}(
    sθ * cϕ,
    sθ * sϕ,
    cθ,
    r * cθ * cϕ,
    r * cθ * sϕ,
    -r * sθ,
    -r * sθ * sϕ,
    r * sθ * cϕ,
    zero(T),
  )
end

@inline function _face_area_vector_from_jacobian(
  jacobian::SMatrix{2,2,T,4}, axis::Int
) where {T}
  if axis == 1
    return SVector{2,T}(jacobian[2, 2], -jacobian[1, 2])
  elseif axis == 2
    return SVector{2,T}(-jacobian[2, 1], jacobian[1, 1])
  end
  throw(ArgumentError("Invalid 2D face axis: $axis"))
end

@inline function _face_area_vector_from_jacobian(
  jacobian::SMatrix{3,3,T,9}, axis::Int
) where {T}
  aξ = SVector{3,T}(jacobian[1, 1], jacobian[2, 1], jacobian[3, 1])
  aη = SVector{3,T}(jacobian[1, 2], jacobian[2, 2], jacobian[3, 2])
  aζ = SVector{3,T}(jacobian[1, 3], jacobian[2, 3], jacobian[3, 3])
  if axis == 1
    return cross(aη, aζ)
  elseif axis == 2
    return cross(aζ, aξ)
  elseif axis == 3
    return cross(aξ, aη)
  end
  throw(ArgumentError("Invalid 3D face axis: $axis"))
end

@inline _face_area_scale(::CoordinateSystemTrait, q::SVector{N,T}) where {N,T} = one(T)
@inline _face_area_scale(::CurvilinearCS, q::SVector{N,T}) where {N,T} = one(T)
@inline _face_area_scale(::CartesianCS, q::SVector{N,T}) where {N,T} = one(T)
@inline _face_area_scale(::SphericalCS, q::SVector{N,T}) where {N,T} = one(T)
@inline _face_area_scale(::CylindricalCS, q::SVector{1,T}) where {T} = T(2π) * abs(q[1])

@inline _face_rotational_radius(::CylindricalCS, q::SVector{2,T}) where {T} = q[1]
@inline _face_rotational_radius(::AxisymmetricCS{:x}, q::SVector{2,T}) where {T} = q[2]
@inline _face_rotational_radius(::AxisymmetricCS{:y}, q::SVector{2,T}) where {T} = q[1]
function _face_rotational_radius(::AxisymmetricCS{Axis}, ::SVector{2,T}) where {Axis,T}
  throw(
    ArgumentError(
      "Unsupported axisymmetric axis `:$Axis`. Supported values are `:x` and `:y`."
    ),
  )
end

@inline function _face_area_scale(
  cs::Union{CylindricalCS,AxisymmetricCS}, q::SVector{2,T}
) where {T}
  return T(2π) * abs(_face_rotational_radius(cs, q))
end

@inline conservation_cell_metric_scale(cs::CoordinateSystemTrait, bt::BasisTrait, q::SVector{N,T}) where {N,T} =
  one(T)
@inline conservation_cell_metric_scale(::CylindricalCS, ::CartesianBasis, q::SVector{1,T}) where {T} =
  T(2π) * abs(q[1])
@inline conservation_cell_metric_scale(::CylindricalCS, ::CartesianBasis, q::SVector{2,T}) where {T} =
  T(2π) * abs(q[1])
@inline conservation_cell_metric_scale(
  ::AxisymmetricCS{:x}, ::CartesianBasis, q::SVector{2,T}
) where {T} = T(2π) * abs(q[2])
@inline conservation_cell_metric_scale(
  ::AxisymmetricCS{:y}, ::CartesianBasis, q::SVector{2,T}
) where {T} = T(2π) * abs(q[1])
@inline conservation_cell_metric_scale(::SphericalCS, ::SphericalBasis, q::SVector{1,T}) where {T} =
  T(4π) * abs(q[1])^2
@inline conservation_cell_metric_scale(::SphericalCS, ::SphericalBasis, q::SVector{2,T}) where {T} =
  T(2π) * abs(q[1])^2 * sin(q[2])
@inline conservation_cell_metric_scale(::SphericalCS, ::SphericalBasis, q::SVector{3,T}) where {T} =
  abs(q[1])^2 * sin(q[2])

@inline function conservation_face_metric_component_scale(cs, bt, q::SVector{N,T}) where {N,T}
  scale = conservation_cell_metric_scale(cs, bt, q)
  return SVector{N,T}(ntuple(_ -> scale, Val(N)))
end

@inline function conservation_face_metric_component_scale(::SphericalCS, ::SphericalBasis, q::SVector{1,T}) where {T}
  r = abs(q[1])
  return SVector{1,T}(T(4π) * r^2)
end

@inline function conservation_face_metric_component_scale(::SphericalCS, ::SphericalBasis, q::SVector{2,T}) where {T}
  r = abs(q[1])
  θ = q[2]
  sθ = sin(θ)
  return SVector{2,T}(T(2π) * r^2 * sθ, T(2π) * r * sθ)
end

@inline function conservation_face_metric_component_scale(::SphericalCS, ::SphericalBasis, q::SVector{3,T}) where {T}
  r = abs(q[1])
  θ = q[2]
  sθ = sin(θ)
  return SVector{3,T}(r^2 * sθ, r * sθ, r)
end

@inline _face_outward_sign(side::Symbol, ::Type{T}) where {T} =
  side === :hi ? one(T) : -one(T)

@inline function _face_geometry(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  if N != 2 && N != 3
    throw(
      ArgumentError(
        "`face_area` and `outward_face_normal` are only supported for 2D and 3D unified grids.",
      ),
    )
  end
  axis, side = _face_loc_axis_side(Val(N), loc)
  Icell = CartesianIndex(idx)
  q = face_coordinate(grid, idx, loc)

  Fq = _face_forward_jacobian(grid, Icell, axis, side)

  cs = coordinate_system(grid)
  A = _face_coord_to_cartesian_jacobian(cs, q)
  Jx = SMatrix{N,N,T,N * N}(Tuple(A * Fq))
  area_vec_plus = _face_area_vector_from_jacobian(Jx, axis)
  outward_vec = area_vec_plus * _face_outward_sign(side, T)
  base_area = norm(outward_vec)
  normal = base_area > zero(T) ? outward_vec / base_area : zero(SVector{N,T})
  area = base_area * _face_area_scale(cs, q)
  x = _face_coord_to_cartesian(cs, q)

  return (; normal=normal, area=area, cartesian_coordinate=x, mapped_coordinate=q)
end

@inline function _face_flux_cache_index(
  idx::NTuple{N,Int}, axis::Int, side::Symbol
) where {N}
  return ntuple(d -> d == axis ? idx[d] + (side === :hi ? 0 : -1) : idx[d], N)
end

"""
    face_flux_geometry(grid::Union{MappedGrid,DiscreteGrid}, idx, loc)

Return the solver-facing face geometry for the face selected by `loc` at cell
index `idx`.

The returned metric vector and normal are expressed in the grid's physical
basis, not in Cartesian embedding coordinates.

For a face on computational axis `alpha`, the returned `metric_vector` is the
outward-oriented physical conserved face metric vector `S^alpha_phys` consumed
by inviscid flux assembly. Coordinate-system measure factors such as
cylindrical/axisymmetric `2πr` and spherical basis factors are included here.

`face_flux_geometry` is intentionally unavailable for `OrthogonalGrid`; those
paths should use [`cellvolume`](@ref), [`face_area`](@ref),
[`face_coordinate`](@ref), and [`basis_transfer_matrix`](@ref) directly.
"""
@inline function face_flux_geometry(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex, loc::Symbol
)
  face_flux_geometry(grid, idx.I, loc)
end

@inline function face_flux_geometry(
  grid::Union{MappedGrid{N},DiscreteGrid{N}}, idx::NTuple{N,Int}, loc::Symbol
) where {N}
  axis, side = _face_loc_axis_side(Val(N), loc)
  face_idx = _face_flux_cache_index(idx, axis, side)
  q = face_coordinate(grid, idx, loc)
  Ghat = face_metrics(grid)[axis].conserved[face_idx...].jacobian_matrix
  T = promote_type(eltype(q), eltype(Ghat))
  metric_row = SVector{N,T}(ntuple(j -> Ghat[axis, j], N))
  qvec = SVector{N,T}(q)
  component_scale = conservation_face_metric_component_scale(coordinate_system(grid), basis_trait(grid), qvec)
  metric_vector = (side === :hi ? one(T) : -one(T)) * (component_scale .* metric_row)
  area = norm(metric_vector)
  normal = area > zero(T) ? metric_vector / area : zero(metric_vector)
  return FaceFluxGeometry{N,T}(qvec, metric_vector, area, normal)
end

function face_flux_geometry(::AbstractOrthogonalGrid, idx, loc::Symbol)
  _ = idx
  _ = loc
  throw(
    ArgumentError(
      "`face_flux_geometry` is undefined for `AbstractOrthogonalGrid`; use `cellvolume`, `face_area`, `face_coordinate`, and `basis_transfer_matrix` instead.",
    ),
  )
end

"""
    face_area(grid::Union{MappedGrid,DiscreteGrid}, idx, loc::Symbol)

Return the physical area (3D) or boundary measure (2D) of the face at cell index
`idx` and face selector `loc`.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: Cell index as `CartesianIndex` or `NTuple{N,Int}`.
  - `loc`: Face selector symbol (for example `:ilo`, `:ihi`, `:jlo`, `:jhi`, `:klo`, `:khi`).

"""
@inline function face_area(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex, loc::Symbol
)
  face_area(grid, idx.I, loc)
end
@inline function face_area(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  geom = _face_geometry(grid, idx, loc)
  return geom.area
end

"""
    outward_face_normal(grid::Union{MappedGrid,DiscreteGrid}, idx, loc::Symbol)

Return the outward Cartesian unit normal vector for the face at cell index `idx`
and face selector `loc`.

This is the embedded boundary or surface normal. It is distinct from
`face_flux_geometry(...).normal`, which is the solver-facing flux normal in the
grid's physical basis.

# Arguments
  - `grid`: Mapped or discrete unified grid.
  - `idx`: Cell index as `CartesianIndex` or `NTuple{N,Int}`.
  - `loc`: Face selector symbol (for example `:ilo`, `:ihi`, `:jlo`, `:jhi`, `:klo`, `:khi`).

"""
@inline function outward_face_normal(
  grid::Union{MappedGrid,DiscreteGrid}, idx::CartesianIndex, loc::Symbol
)
  outward_face_normal(grid, idx.I, loc)
end
@inline function outward_face_normal(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}, loc::Symbol
) where {N,T}
  geom = _face_geometry(grid, idx, loc)
  return geom.normal
end

@inline function face_area(grid::AbstractOrthogonalGrid, dim::Int, idx::CartesianIndex)
  face_area(grid, dim, idx.I)
end
@inline function face_area(grid::OrthogonalGrid{N}, dim::Int, idx::NTuple{N,Int}) where {N}
  1 <= dim <= N ||
    throw(ArgumentError("face axis dim=$dim is invalid for $N-D orthogonal grid"))
  grid.face_areas[dim][idx...]
end

@inline function face_area(grid::AbstractOrthogonalGrid, idx::CartesianIndex, loc::Symbol)
  face_area(grid, idx.I, loc)
end
@inline function face_area(
  grid::OrthogonalGrid{N}, idx::NTuple{N,Int}, loc::Symbol
) where {N}
  axis, side = _face_loc_axis_side(Val(N), loc)
  ioff = side === :hi ? 1 : 0
  i = ntuple(d -> (d == axis ? idx[d] + ioff : idx[d]), N)
  face_area(grid, axis, i)
end

function outward_face_normal(::AbstractOrthogonalGrid, idx, loc::Symbol)
  throw(ArgumentError("`outward_face_normal` is undefined for `AbstractOrthogonalGrid`."))
end

@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::NTuple{N,Int}
) where {N,T}
  return abs(_cell_forward_metric_at(grid, idx).J)
end
@inline function _jacobian_volume_factor(
  grid::Union{MappedGrid{N,T},DiscreteGrid{N,T}}, idx::Tuple{Vararg{Real,N}}
) where {N,T}
  return abs(T(det(_continuous_forward_jacobian(grid, idx))))
end

@inline _radial_centroid_1d(
  grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::NTuple{1,Int}
) = grid.centroid_coordinates[1][idx...]
@inline _radial_centroid_1d(
  grid::Union{MappedGrid{1},DiscreteGrid{1}}, idx::Tuple{Vararg{Real,1}}
) = _continuous_coord(grid, idx)[1]

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
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::NTuple{2,Int},
)
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::CylindricalCS,
  ::CartesianBasis,
  grid::Union{MappedGrid{2},DiscreteGrid{2}},
  idx::Tuple{Vararg{Real,2}},
)
  return conservation_cell_metric_scale(CylindricalCS(), CartesianBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::NTuple{2,Int},
) where {Axis}
  return conservation_cell_metric_scale(cs, CartesianBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  cs::AxisymmetricCS{Axis},
  ::CartesianBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
) where {Axis}
  return conservation_cell_metric_scale(cs, CartesianBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{1,Int}
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,1}},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::NTuple{3,Int},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

@inline function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::Union{MappedGrid{3},DiscreteGrid{3}},
  idx::Tuple{Vararg{Real,3}},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  ::SphericalCS, ::SphericalBasis, grid::AbstractMappedOrDiscreteGrid, idx::NTuple{2,Int}
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), centroid(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
end

function _cellvolume_dispatch(
  ::SphericalCS,
  ::SphericalBasis,
  grid::AbstractMappedOrDiscreteGrid,
  idx::Tuple{Vararg{Real,2}},
)
  return conservation_cell_metric_scale(SphericalCS(), SphericalBasis(), _continuous_coord(grid, idx)) *
         _jacobian_volume_factor(grid, idx)
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

function cellvolumes(grid::Union{MappedGrid{N},DiscreteGrid{N}}; include_halo::Bool=false) where {N}
  volumes = KernelAbstractions.zeros(
    grid.backend, eltype(grid), size(grid.iterators.cell.full); unified=false
  )
  _compute_unified_cell_volumes!(volumes, grid, Val(N))
  include_halo && return volumes
  @views return volumes[grid.iterators.cell.domain]
end

function cellvolumes(grid::OrthogonalGrid; include_halo::Bool=false)
  include_halo && return grid.cell_volumes
  @views return grid.cell_volumes[grid.iterators.cell.domain]
end

function _compute_unified_cell_volumes!(volumes, grid, ::Val{N}) where {N}
  backend = KernelAbstractions.get_backend(volumes)
  kernel = _compute_unified_cell_volumes_kernel!(backend)
  metrics = cell_metrics(grid)
  kernel(
    volumes,
    metrics.forward,
    grid.centroid_coordinates,
    coordinate_system(grid),
    basis_trait(grid),
    grid.iterators.cell.full,
    Val(N);
    ndrange=length(grid.iterators.cell.full),
  )
  KernelAbstractions.synchronize(backend)
  return volumes
end

@kernel function _compute_unified_cell_volumes_kernel!(
  volumes,
  cell_forward_metrics,
  centroid_coordinates,
  coordinate_system,
  basis,
  local_domain,
  ::Val{N},
) where {N}
  idx = @index(Global, Linear)
  I = local_domain[idx]
  q = SVector{N}(ntuple(d -> centroid_coordinates[d][I], N))
  scale = conservation_cell_metric_scale(coordinate_system, basis, q)
  volumes[I] = scale * abs(cell_forward_metrics[I].J)
end

cellsize(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.domain)
cellsize(grid::OrthogonalGrid) = size(grid.iterators.cell.domain)

cellsize_withhalo(grid::Union{MappedGrid,DiscreteGrid}) = size(grid.iterators.cell.full)
cellsize_withhalo(grid::OrthogonalGrid) = size(grid.iterators.cell.full)

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
