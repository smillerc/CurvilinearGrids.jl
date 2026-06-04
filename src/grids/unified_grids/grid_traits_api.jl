#
# Grid traits API
#

#
# Coordinate and basis traits
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
#
# Basis transfer
#

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
#
# Element type and compatibility constructors
#

Base.eltype(::MappedGrid{N,T}) where {N,T} = T
Base.eltype(::DiscreteGrid{N,T}) where {N,T} = T

# Backward-compatible identity constructor for orthogonal grids.
OrthogonalGrid(grid::OrthogonalGrid) = grid
