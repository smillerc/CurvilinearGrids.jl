#
# Face geometry
#

#
# Face flux geometry type
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
#
# Face location helpers
#

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
#
# Face coordinate accessors
#

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
#
# Face embedding transforms
#

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
#
# Coordinate-system measure scales
#

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
#
# Face geometry accessors
#

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
