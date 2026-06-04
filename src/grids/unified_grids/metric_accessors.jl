#
# Metric accessors
#

#
# Cell metric accessors
#

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
#
# Face metric coefficients
#

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
#
# Jacobian accessors
#

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
