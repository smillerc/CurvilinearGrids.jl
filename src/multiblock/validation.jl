"""
    validate_multiblock!(mb::MultiBlockMesh)

Validate multi-block topology and geometric abutting constraints.

# Arguments
  - `mb`: Multi-block mesh container.

# Returns
`mb` on success.
"""
function validate_multiblock!(mb::MultiBlockMesh{N}) where {N}
  _validate_interface_faces_unique!(mb, Val(N))
  for (iface_id, iface) in pairs(mb.interfaces)
    _validate_interface!(mb, iface_id, iface, Val(N))
  end
  return mb
end

function _validate_interface_faces_unique!(mb::MultiBlockMesh{N}, ::Val{N}) where {N}
  nblocks = length(mb.blocks)
  claims = Dict{Tuple{Int,Int,Symbol},Int}()

  for (iface_id, iface) in pairs(mb.interfaces)
    for face in (iface.left, iface.right)
      _validate_face_reference(face, nblocks, N, "interface[$iface_id]")
      key = _face_key(face)
      if haskey(claims, key)
        throw(
          ArgumentError(
            "Face $key is multiply assigned across interfaces (existing=interface[$(claims[key])], new=interface[$iface_id]).",
          ),
        )
      end
      claims[key] = iface_id
    end
  end
  return nothing
end

function _validate_face_reference(
  face::BlockFace{N}, nblocks::Int, ndims::Int, where_str
) where {N}
  if face.block_id < 1 || face.block_id > nblocks
    throw(
      ArgumentError(
        "Invalid block id $(face.block_id) in $where_str. Expected in 1:$nblocks."
      ),
    )
  end
  if face.axis < 1 || face.axis > ndims
    throw(ArgumentError("Invalid axis $(face.axis) in $where_str for N=$ndims."))
  end
  if face.side !== :min && face.side !== :max
    throw(ArgumentError("Invalid side $(face.side) in $where_str."))
  end
  return nothing
end

function _validate_interface!(
  mb::MultiBlockMesh{N}, iface_id::Int, iface::BlockInterface{N}, ::Val{N}
) where {N}
  left_block = mb.blocks[iface.left.block_id]
  right_block = mb.blocks[iface.right.block_id]

  left_tangential = _face_tangential_ranges(left_block, iface.left, Val(N))
  right_tangential = _face_tangential_ranges(right_block, iface.right, Val(N))

  left_sizes = ntuple(i -> length(left_tangential[i]), N - 1)
  right_sizes = ntuple(i -> length(right_tangential[i]), N - 1)

  for i in 1:(N - 1)
    right_dim = iface.permutation[i]
    if left_sizes[i] != right_sizes[right_dim]
      throw(
        ArgumentError(
          "Interface[$iface_id] tangential size mismatch: left dim $i has $(left_sizes[i]), right perm dim $right_dim has $(right_sizes[right_dim]).",
        ),
      )
    end
  end

  left_interior, right_interior = _paired_interface_boundary_indices(mb, iface, Val(N))
  _check_abutting_geometry!(mb, iface_id, iface, left_interior, right_interior, Val(N))
  return nothing
end

function _check_abutting_geometry!(
  mb::MultiBlockMesh{N,T},
  iface_id::Int,
  iface::BlockInterface{N},
  left_interior::Vector{CartesianIndex{N}},
  right_interior::Vector{CartesianIndex{N}},
  ::Val{N},
) where {N,T}
  left_block = mb.blocks[iface.left.block_id]
  right_block = mb.blocks[iface.right.block_id]

  left_points = Vector{SVector{N,T}}(undef, length(left_interior))
  right_points = Vector{SVector{N,T}}(undef, length(right_interior))

  max_distance = zero(T)
  for i in eachindex(left_interior)
    ql = _face_point_native(left_block, left_interior[i], iface.left, T, Val(N))
    qr = _face_point_native(right_block, right_interior[i], iface.right, T, Val(N))
    xl = _coord_to_cartesian(coordinate_system(left_block), ql)
    xr = _coord_to_cartesian(coordinate_system(right_block), qr)
    left_points[i] = xl
    right_points[i] = xr
    max_distance = max(max_distance, norm(xl - xr))
  end

  _check_aabb_overlap!(iface_id, left_points, right_points, iface.tolerance)
  if max_distance > iface.tolerance
    throw(
      ArgumentError(
        "Interface[$iface_id] is not abutting: max point distance $max_distance exceeds tolerance $(iface.tolerance).",
      ),
    )
  end

  normal_tol = max(sqrt(iface.tolerance), eps(T))
  for i in eachindex(left_interior)
    nl = _face_outward_normal_cartesian(left_block, left_interior[i], iface.left, Val(N))
    nr = _face_outward_normal_cartesian(right_block, right_interior[i], iface.right, Val(N))
    if dot(nl, nr) > -one(T) + normal_tol
      throw(
        ArgumentError(
          "Interface[$iface_id] normals are not opposing at sample $i (dot=$(dot(nl, nr)))."
        ),
      )
    end
  end
  return nothing
end

function _check_aabb_overlap!(iface_id::Int, left_points, right_points, tol::Real)
  N = length(first(left_points))
  left_min = fill(typemax(eltype(first(left_points))), N)
  left_max = fill(typemin(eltype(first(left_points))), N)
  right_min = fill(typemax(eltype(first(right_points))), N)
  right_max = fill(typemin(eltype(first(right_points))), N)

  for p in left_points
    for d in 1:N
      left_min[d] = min(left_min[d], p[d])
      left_max[d] = max(left_max[d], p[d])
    end
  end
  for p in right_points
    for d in 1:N
      right_min[d] = min(right_min[d], p[d])
      right_max[d] = max(right_max[d], p[d])
    end
  end

  for d in 1:N
    if left_max[d] < right_min[d] - tol || right_max[d] < left_min[d] - tol
      throw(
        ArgumentError("Interface[$iface_id] bounding boxes do not abut along dimension $d.")
      )
    end
  end
  return nothing
end

function _face_tangential_ranges(
  grid::AbstractMappedOrDiscreteGrid, face::BlockFace{N}, ::Val{N}
) where {N}
  domain_ranges = Tuple(grid.iterators.cell.domain.indices)
  return ntuple(i -> domain_ranges[i < face.axis ? i : i + 1], N - 1)
end

function _face_tangential_ranges(grid::OrthogonalGrid{N}, face::BlockFace{N}, ::Val{N}) where {N}
  domain_ranges = Tuple(grid.iterators.cell.domain.indices)
  return ntuple(i -> domain_ranges[i < face.axis ? i : i + 1], N - 1)
end

function _paired_interface_boundary_indices(
  mb::MultiBlockMesh{N}, iface::BlockInterface{N}, ::Val{N}
) where {N}
  left_block = mb.blocks[iface.left.block_id]
  right_block = mb.blocks[iface.right.block_id]

  left_tangential = _face_tangential_ranges(left_block, iface.left, Val(N))
  right_tangential = _face_tangential_ranges(right_block, iface.right, Val(N))

  left_axis_idx = _face_axis_index(
    Tuple(left_block.iterators.cell.domain.indices), iface.left, 0
  )
  right_axis_idx = _face_axis_index(
    Tuple(right_block.iterators.cell.domain.indices), iface.right, 0
  )

  left_interior = CartesianIndex{N}[]
  right_interior = CartesianIndex{N}[]

  if N == 1
    push!(left_interior, CartesianIndex(left_axis_idx))
    push!(right_interior, CartesianIndex(right_axis_idx))
    return left_interior, right_interior
  end

  for left_t in CartesianIndices(left_tangential)
    left_tvals = left_t.I
    right_tvals = _map_tangential_indices(
      left_tvals, left_tangential, right_tangential, iface.permutation, iface.flips, Val(N)
    )
    push!(
      left_interior,
      _build_face_cartesian_index(left_tvals, iface.left.axis, left_axis_idx, Val(N)),
    )
    push!(
      right_interior,
      _build_face_cartesian_index(right_tvals, iface.right.axis, right_axis_idx, Val(N)),
    )
  end
  return left_interior, right_interior
end

@inline function _face_axis_index(
  domain_ranges::NTuple{N,UnitRange{Int}}, face::BlockFace{N}, ghost_layer::Int
) where {N}
  axis_range = domain_ranges[face.axis]
  if face.side === :min
    return first(axis_range) - ghost_layer
  end
  return last(axis_range) + ghost_layer
end

@inline function _build_face_cartesian_index(
  tangential::NTuple{M,Int}, axis::Int, axis_index::Int, ::Val{N}
) where {M,N}
  if M != N - 1
    throw(
      ArgumentError(
        "Invalid tangential tuple length $M for face in dimension $N (expected $(N - 1))."
      ),
    )
  end
  full = ntuple(i -> begin
    if i == axis
      return axis_index
    end
    ti = i < axis ? i : i - 1
    return tangential[ti]
  end, N)
  return CartesianIndex(full)
end

@inline function _map_tangential_indices(
  left_tvals::Tuple,
  left_ranges::Tuple,
  right_ranges::Tuple,
  permutation::Tuple,
  flips::Tuple,
  ::Val{N},
) where {N}
  n_tangential = N - 1
  if length(left_tvals) != n_tangential ||
    length(left_ranges) != n_tangential ||
    length(right_ranges) != n_tangential ||
    length(permutation) != n_tangential ||
    length(flips) != n_tangential
    throw(ArgumentError("Tangential tuple lengths are inconsistent for N=$N."))
  end
  right_values = Vector{Int}(undef, n_tangential)
  for i in 1:(N - 1)
    right_dim = permutation[i]
    offset = left_tvals[i] - first(left_ranges[i])
    if flips[i]
      offset = (length(right_ranges[right_dim]) - 1) - offset
    end
    right_values[right_dim] = first(right_ranges[right_dim]) + offset
  end
  return Tuple(right_values)
end

@inline function _as_coord_svector(x, ::Type{T}, N::Int) where {T}
  return SVector{N,T}(ntuple(i -> T(x[i]), N))
end

@inline function _coord_to_cartesian(::CartesianCS, q::SVector{N,T}) where {N,T}
  return q
end
@inline function _coord_to_cartesian(::CurvilinearCS, q::SVector{N,T}) where {N,T}
  return q
end
@inline function _coord_to_cartesian(::CoordinateSystemTrait, q::SVector{N,T}) where {N,T}
  return q
end
@inline function _coord_to_cartesian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return @SVector [r * sθ * cϕ, r * sθ * sϕ, r * cθ]
end
@inline function _coord_to_cartesian(::SphericalCS, q::SVector{2,T}) where {T}
  r, θ = q
  return @SVector [r * sin(θ), r * cos(θ)]
end

@inline function _coord_to_cartesian_jacobian(::CartesianCS, q::SVector{N,T}) where {N,T}
  return one(SMatrix{N,N,T})
end
@inline function _coord_to_cartesian_jacobian(::CurvilinearCS, q::SVector{N,T}) where {N,T}
  return one(SMatrix{N,N,T})
end
@inline function _coord_to_cartesian_jacobian(
  ::CoordinateSystemTrait, q::SVector{N,T}
) where {N,T}
  return one(SMatrix{N,N,T})
end
@inline function _coord_to_cartesian_jacobian(::SphericalCS, q::SVector{3,T}) where {T}
  r, θ, ϕ = q
  sθ = sin(θ)
  cθ = cos(θ)
  sϕ = sin(ϕ)
  cϕ = cos(ϕ)
  return @SMatrix [
    sθ*cϕ r*cθ*cϕ -r*sθ*sϕ
    sθ*sϕ r*cθ*sϕ r*sθ*cϕ
    cθ -r*sθ 0
  ]
end
@inline function _coord_to_cartesian_jacobian(::SphericalCS, q::SVector{2,T}) where {T}
  r, θ = q
  sθ = sin(θ)
  cθ = cos(θ)
  return @SMatrix [
    sθ r * cθ
    cθ -r * sθ
  ]
end

@inline function _face_point_native(
  grid::AbstractMappedOrDiscreteGrid, idx::CartesianIndex{N}, face::BlockFace{N}, ::Type{T}, ::Val{N}
) where {N,T}
  q = face_coordinate(grid, idx.I, _block_face_location(face, Val(N)))
  return _as_coord_svector(q, T, N)
end

@inline function _face_point_native(
  grid::AbstractOrthogonalGrid, idx::CartesianIndex{N}, face::BlockFace{N}, ::Type{T}, ::Val{N}
) where {N,T}
  q = face_coordinate(grid, idx.I, _block_face_location(face, Val(N)))
  return _as_coord_svector(q, T, N)
end

@inline function _face_outward_normal_cartesian(
  grid::AbstractMappedOrDiscreteGrid, idx::CartesianIndex{N}, face::BlockFace{N}, ::Val{N}
) where {N}
  T = eltype(grid)
  ξ = _face_computational_coordinate(idx, face, T, Val(N))
  q = _as_coord_svector(coord(grid, ξ), T, N)
  Jq = forward_cell_metrics(grid, ξ).jacobian_matrix
  A = _coord_to_cartesian_jacobian(coordinate_system(grid), q)
  Jx = A * Jq
  G = inv(Jx)
  n = SVector{N,T}(ntuple(i -> G[face.axis, i], N))
  n_out = face.side === :min ? -n : n
  n_norm = norm(n_out)
  if n_norm <= eps(T)
    throw(ArgumentError("Degenerate interface normal encountered at index $idx."))
  end
  return n_out / n_norm
end

@inline function _face_outward_normal_cartesian(
  grid::AbstractOrthogonalGrid, idx::CartesianIndex{N}, face::BlockFace{N}, ::Val{N}
) where {N}
  T = eltype(grid)
  q = _face_point_native(grid, idx, face, T, Val(N))
  A = _coord_to_cartesian_jacobian(coordinate_system(grid), q)
  n = SVector{N,T}(ntuple(i -> A[i, face.axis], N))
  n_out = face.side === :min ? -n : n
  n_norm = norm(n_out)
  if n_norm <= eps(T)
    throw(ArgumentError("Degenerate interface normal encountered at index $idx."))
  end
  return n_out / n_norm
end

@inline function _block_face_location(face::BlockFace{N}, ::Val{N}) where {N}
  if face.axis == 1
    return face.side === :min ? :ilo : :ihi
  elseif face.axis == 2 && N >= 2
    return face.side === :min ? :jlo : :jhi
  elseif face.axis == 3 && N >= 3
    return face.side === :min ? :klo : :khi
  end
  throw(ArgumentError("Invalid block face `(axis=$(face.axis), side=$(face.side))` for N=$N."))
end

@inline function _face_computational_coordinate(
  idx::CartesianIndex{N}, face::BlockFace{N}, ::Type{T}, ::Val{N}
) where {N,T}
  return ntuple(
    d -> begin
      if d == face.axis
        return face.side === :max ? T(idx.I[d]) + T(0.5) : T(idx.I[d]) - T(0.5)
      end
      return T(idx.I[d])
    end, N
  )
end
