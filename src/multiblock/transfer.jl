"""
    exchange_interface!(mb, iface_id, fields; field_kind=:scalar, direction=:both, depth=1)

Exchange one interface for cell-centered fields.

# Arguments
  - `mb`: Multi-block mesh.
  - `iface_id`: Interface index.
  - `fields`: Block-aligned vector of arrays (`fields[block_id]`).

# Keywords
  - `field_kind`: `:scalar`, `:vector`, or `:tensor`.
  - `direction`: `:left_to_right`, `:right_to_left`, or `:both`.
  - `depth`: Number of halo layers to fill from the opposing block interior.

# Returns
`fields`.
"""
function exchange_interface!(
  mb::MultiBlockMesh{N,T},
  iface_id::Int,
  fields;
  field_kind::Symbol=:scalar,
  direction::Symbol=:both,
  depth::Integer=1,
) where {N,T}
  if iface_id < 1 || iface_id > length(mb.interfaces)
    throw(ArgumentError("Invalid `iface_id=$iface_id`."))
  end
  if field_kind ∉ (:scalar, :vector, :tensor)
    throw(ArgumentError("Invalid `field_kind=$field_kind`."))
  end
  if direction ∉ (:left_to_right, :right_to_left, :both)
    throw(ArgumentError("Invalid `direction=$direction`."))
  end
  depth >= 1 || throw(ArgumentError("`depth` must be positive."))

  _validate_fields_container(mb, fields)
  _validate_exchange_depth(mb, depth)
  iface = mb.interfaces[iface_id]

  if depth > 1
    for layer in 1:depth
      _exchange_interface_layer!(mb, iface, fields, field_kind, direction, layer, Val(N), T)
    end
    return fields
  end

  cache = _ensure_interface_cache!(mb, iface_id)

  left_values = fields[iface.left.block_id]
  right_values = fields[iface.right.block_id]

  if direction === :left_to_right || direction === :both
    _exchange_one_direction!(
      left_values,
      right_values,
      cache.left_interior,
      cache.right_ghost,
      cache.left_to_right_transform,
      field_kind,
      Val(N),
      T,
    )
  end
  if direction === :right_to_left || direction === :both
    _exchange_one_direction!(
      right_values,
      left_values,
      cache.right_interior,
      cache.left_ghost,
      cache.right_to_left_transform,
      field_kind,
      Val(N),
      T,
    )
  end
  return fields
end

"""
    exchange_all_interfaces!(mb, fields; field_kind=:scalar, direction=:both, depth=1)

Exchange all interfaces for cell-centered fields.

# Arguments
  - `mb`: Multi-block mesh.
  - `fields`: Block-aligned vector of arrays (`fields[block_id]`).

# Keywords
  - `field_kind`: `:scalar`, `:vector`, or `:tensor`.
  - `direction`: `:left_to_right`, `:right_to_left`, or `:both`.
  - `depth`: Number of halo layers to fill from the opposing block interior.

# Returns
`fields`.
"""
function exchange_all_interfaces!(
  mb::MultiBlockMesh,
  fields;
  field_kind::Symbol=:scalar,
  direction::Symbol=:both,
  depth::Integer=1,
)
  for iface_id in eachindex(mb.interfaces)
    exchange_interface!(
      mb, iface_id, fields; field_kind=field_kind, direction=direction, depth=depth
    )
  end
  return fields
end

function _validate_fields_container(mb::MultiBlockMesh, fields)
  if length(fields) != length(mb.blocks)
    throw(
      ArgumentError(
        "Expected one field array per block (expected=$(length(mb.blocks)), got=$(length(fields))).",
      ),
    )
  end
  return nothing
end

function _validate_exchange_depth(mb::MultiBlockMesh, depth::Integer)
  for (block_id, block) in pairs(mb.blocks)
    if depth > block.nhalo
      throw(
        ArgumentError(
          "`depth=$depth` exceeds halo width $(block.nhalo) for block $block_id.",
        ),
      )
    end
  end
  return nothing
end

function _exchange_interface_layer!(
  mb::MultiBlockMesh{N},
  iface::BlockInterface{N},
  fields,
  field_kind::Symbol,
  direction::Symbol,
  layer::Int,
  ::Val{N},
  ::Type{T},
) where {N,T}
  left_values = fields[iface.left.block_id]
  right_values = fields[iface.right.block_id]
  left_interior, right_interior, left_ghost, right_ghost =
    _paired_interface_layer_indices(mb, iface, layer, Val(N))
  left_to_right_transform, right_to_left_transform = if field_kind === :scalar
    nothing, nothing
  else
    _interface_layer_transforms(mb, iface, left_interior, right_interior, Val(N), T)
  end

  if direction === :left_to_right || direction === :both
    _exchange_one_direction!(
      left_values,
      right_values,
      left_interior,
      right_ghost,
      left_to_right_transform,
      field_kind,
      Val(N),
      T,
    )
  end
  if direction === :right_to_left || direction === :both
    _exchange_one_direction!(
      right_values,
      left_values,
      right_interior,
      left_ghost,
      right_to_left_transform,
      field_kind,
      Val(N),
      T,
    )
  end
  return nothing
end

function _paired_interface_layer_indices(
  mb::MultiBlockMesh{N}, iface::BlockInterface{N}, layer::Int, ::Val{N}
) where {N}
  left_block = mb.blocks[iface.left.block_id]
  right_block = mb.blocks[iface.right.block_id]

  left_domain = Tuple(left_block.iterators.cell.domain.indices)
  right_domain = Tuple(right_block.iterators.cell.domain.indices)
  left_tangential = _face_tangential_ranges(left_block, iface.left, Val(N))
  right_tangential = _face_tangential_ranges(right_block, iface.right, Val(N))

  source_layer = -(layer - 1)
  left_axis_idx = _face_axis_index(left_domain, iface.left, source_layer)
  right_axis_idx = _face_axis_index(right_domain, iface.right, source_layer)
  left_ghost_axis = _face_axis_index(left_domain, iface.left, layer)
  right_ghost_axis = _face_axis_index(right_domain, iface.right, layer)

  left_interior = CartesianIndex{N}[]
  right_interior = CartesianIndex{N}[]
  left_ghost = CartesianIndex{N}[]
  right_ghost = CartesianIndex{N}[]

  if N == 1
    push!(left_interior, CartesianIndex(left_axis_idx))
    push!(right_interior, CartesianIndex(right_axis_idx))
    push!(left_ghost, CartesianIndex(left_ghost_axis))
    push!(right_ghost, CartesianIndex(right_ghost_axis))
    return left_interior, right_interior, left_ghost, right_ghost
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
    push!(
      left_ghost,
      _build_face_cartesian_index(left_tvals, iface.left.axis, left_ghost_axis, Val(N)),
    )
    push!(
      right_ghost,
      _build_face_cartesian_index(right_tvals, iface.right.axis, right_ghost_axis, Val(N)),
    )
  end
  return left_interior, right_interior, left_ghost, right_ghost
end

function _interface_layer_transforms(
  mb::MultiBlockMesh{N},
  iface::BlockInterface{N},
  left_interior,
  right_interior,
  ::Val{N},
  ::Type{T},
) where {N,T}
  left_block = mb.blocks[iface.left.block_id]
  right_block = mb.blocks[iface.right.block_id]
  left_to_right = Vector{SMatrix{N,N,T}}(undef, length(left_interior))
  right_to_left = Vector{SMatrix{N,N,T}}(undef, length(left_interior))

  for i in eachindex(left_interior)
    ql = _as_coord_svector(centroid(left_block, left_interior[i]), T, N)
    qr = _as_coord_svector(centroid(right_block, right_interior[i]), T, N)
    left_to_right[i] = _basis_transfer_matrix(left_block, ql, right_block, qr, Val(N), T)
    right_to_left[i] = _basis_transfer_matrix(right_block, qr, left_block, ql, Val(N), T)
  end
  return left_to_right, right_to_left
end

function _exchange_one_direction!(
  donor_values,
  receiver_values,
  donor_indices,
  receiver_indices,
  transforms,
  field_kind::Symbol,
  ::Val{N},
  ::Type{T},
) where {N,T}
  for i in eachindex(donor_indices)
    src = donor_indices[i]
    dst = receiver_indices[i]
    if field_kind === :scalar
      receiver_values[dst] = donor_values[src]
    elseif field_kind === :vector
      vsrc = _as_vector_value(donor_values[src], Val(N), T)
      vdst = transforms[i] * vsrc
      receiver_values[dst] = _coerce_like(receiver_values[dst], vdst)
    else
      tsrc = _as_tensor_value(donor_values[src], Val(N), T)
      A = transforms[i]
      tdst = A * tsrc * transpose(A)
      receiver_values[dst] = _coerce_like(receiver_values[dst], tdst)
    end
  end
  return nothing
end

@inline function _as_vector_value(x, ::Val{N}, ::Type{T}) where {N,T}
  return SVector{N,T}(ntuple(i -> T(x[i]), N))
end

@inline function _as_tensor_value(x, ::Val{N}, ::Type{T}) where {N,T}
  return SMatrix{N,N,T}(ntuple(i -> T(x[i]), N * N))
end

@inline _coerce_like(::Any, value) = value
@inline function _coerce_like(example::SVector{N,T}, value::SVector{N,T}) where {N,T}
  return value
end
@inline function _coerce_like(
  example::SMatrix{N,N,T,NN}, value::SMatrix{N,N,T,NN}
) where {N,T,NN}
  return value
end
@inline function _coerce_like(example::NTuple{N,T}, value::SVector{N,T}) where {N,T}
  return Tuple(value)
end

@inline function _extract_tangential_coordinate(
  ξ::NTuple{N,T}, axis::Int, ::Val{N}
) where {N,T}
  return ntuple(i -> ξ[i < axis ? i : i + 1], N - 1)
end

@inline function _inverse_interface_map(
  permutation::NTuple{M,Int}, flips::NTuple{M,Bool}
) where {M}
  inv_perm = ntuple(j -> findfirst(==(j), permutation), M)
  inv_flips = ntuple(j -> flips[inv_perm[j]], M)
  return inv_perm, inv_flips
end

@inline function _map_tangential_coordinates(
  source_tangential::Tuple,
  source_ranges::Tuple,
  target_ranges::Tuple,
  permutation::Tuple,
  flips::Tuple,
  ::Val{N},
  ::Type{T},
) where {N,T}
  n_tangential = N - 1
  mapped = Vector{T}(undef, n_tangential)
  for i in 1:n_tangential
    target_dim = permutation[i]
    offset = T(source_tangential[i]) - T(first(source_ranges[i]))
    if flips[i]
      offset = T(length(target_ranges[target_dim]) - 1) - offset
    end
    mapped[target_dim] = T(first(target_ranges[target_dim])) + offset
  end
  return Tuple(mapped)
end

@inline function _face_axis_coordinate(
  domain_ranges::NTuple{N,UnitRange{Int}}, face::BlockFace{N}, ::Type{T}
) where {N,T}
  axis_range = domain_ranges[face.axis]
  if face.side === :min
    return T(first(axis_range)) - T(0.5)
  end
  return T(last(axis_range)) + T(0.5)
end

@inline function _build_face_coordinate(
  tangential::Tuple, axis::Int, axis_coordinate::T, ::Val{N}
) where {N,T}
  return ntuple(i -> begin
    if i == axis
      return axis_coordinate
    end
    ti = i < axis ? i : i - 1
    return T(tangential[ti])
  end, N)
end

"""
    computational_coordinate(mb, iface_id, ξ; from=:left)

Transfer a computational coordinate across one interface.

`ξ` is interpreted in the source block's computational space and mapped to the
opposite block's interface coordinate.

# Keywords
  - `from`: `:left` or `:right`.
"""
function computational_coordinate(
  mb::MultiBlockMesh{N,T}, iface_id::Int, ξ::Tuple{Vararg{Real,N}}; from::Symbol=:left
) where {N,T}
  if iface_id < 1 || iface_id > length(mb.interfaces)
    throw(ArgumentError("Invalid `iface_id=$iface_id`."))
  end
  if from ∉ (:left, :right)
    throw(ArgumentError("Invalid `from=$from`. Expected `:left` or `:right`."))
  end

  iface = mb.interfaces[iface_id]
  if from === :left
    source_face = iface.left
    target_face = iface.right
    permutation = iface.permutation
    flips = iface.flips
  else
    source_face = iface.right
    target_face = iface.left
    permutation, flips = _inverse_interface_map(iface.permutation, iface.flips)
  end

  source_grid = mb.blocks[source_face.block_id]
  target_grid = mb.blocks[target_face.block_id]

  source_ranges = _face_tangential_ranges(source_grid, source_face, Val(N))
  target_ranges = _face_tangential_ranges(target_grid, target_face, Val(N))

  ξsrc = ntuple(i -> T(ξ[i]), N)
  tangential_src = _extract_tangential_coordinate(ξsrc, source_face.axis, Val(N))
  tangential_dst = _map_tangential_coordinates(
    tangential_src, source_ranges, target_ranges, permutation, flips, Val(N), T
  )
  target_axis_coordinate = _face_axis_coordinate(
    Tuple(target_grid.iterators.cell.domain.indices), target_face, T
  )
  return SVector{N,T}(
    _build_face_coordinate(tangential_dst, target_face.axis, target_axis_coordinate, Val(N))
  )
end

@inline function computational_coordinate(
  mb::MultiBlockMesh{N,T}, iface_id::Int, ξ::SVector{N,<:Real}; from::Symbol=:left
) where {N,T}
  return computational_coordinate(mb, iface_id, Tuple(ξ); from=from)
end

"""
    computational_coordinate(mb, iface_id, source_block_id, ξ)

Transfer computational coordinate `ξ` to the block on the opposite side of the
interface selected by `iface_id`.

# Returns
Named tuple `(block_id, coordinate)`.
"""
function computational_coordinate(
  mb::MultiBlockMesh{N,T}, iface_id::Int, source_block_id::Int, ξ::Tuple{Vararg{Real,N}}
) where {N,T}
  if iface_id < 1 || iface_id > length(mb.interfaces)
    throw(ArgumentError("Invalid `iface_id=$iface_id`."))
  end
  iface = mb.interfaces[iface_id]
  if source_block_id == iface.left.block_id
    return (;
      block_id=iface.right.block_id,
      coordinate=computational_coordinate(mb, iface_id, ξ; from=:left),
    )
  elseif source_block_id == iface.right.block_id
    return (;
      block_id=iface.left.block_id,
      coordinate=computational_coordinate(mb, iface_id, ξ; from=:right),
    )
  end
  throw(
    ArgumentError(
      "Block $source_block_id is not part of interface[$iface_id] ($(iface.left.block_id) <-> $(iface.right.block_id)).",
    ),
  )
end

@inline function computational_coordinate(
  mb::MultiBlockMesh{N,T}, iface_id::Int, source_block_id::Int, ξ::SVector{N,<:Real}
) where {N,T}
  return computational_coordinate(mb, iface_id, source_block_id, Tuple(ξ))
end
