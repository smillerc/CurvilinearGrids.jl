"""
    invalidate_interface_caches!(mb::MultiBlockMesh)

Drop all cached interface maps.

# Arguments
  - `mb`: Multi-block mesh.

# Returns
`nothing`.
"""
function invalidate_interface_caches!(mb::MultiBlockMesh)
  fill!(mb.cache.interfaces, nothing)
  return nothing
end

"""
    build_interface_caches!(mb::MultiBlockMesh)

Build (or rebuild) per-interface transfer caches.

# Arguments
  - `mb`: Multi-block mesh.

# Returns
`mb.cache.interfaces`.
"""
function build_interface_caches!(mb::MultiBlockMesh{N,T}) where {N,T}
  for (iface_id, iface) in pairs(mb.interfaces)
    mb.cache.interfaces[iface_id] = _build_interface_cache(mb, iface, Val(N), T)
  end
  return mb.cache.interfaces
end

@inline function _ensure_interface_cache!(
  mb::MultiBlockMesh{N,T}, iface_id::Int
) where {N,T}
  iface = mb.interfaces[iface_id]
  cache = mb.cache.interfaces[iface_id]
  if cache === nothing
    cache = _build_interface_cache(mb, iface, Val(N), T)
    mb.cache.interfaces[iface_id] = cache
    return cache
  end
  left_block = mb.blocks[iface.left.block_id]
  right_block = mb.blocks[iface.right.block_id]
  if cache.left_state_token != _grid_state_token(left_block) ||
    cache.right_state_token != _grid_state_token(right_block)
    cache = _build_interface_cache(mb, iface, Val(N), T)
    mb.cache.interfaces[iface_id] = cache
  end
  return cache
end

function _build_interface_cache(
  mb::MultiBlockMesh{N}, iface::BlockInterface{N}, ::Val{N}, ::Type{T}
) where {N,T}
  left_block = mb.blocks[iface.left.block_id]
  right_block = mb.blocks[iface.right.block_id]

  left_interior, right_interior = _paired_interface_boundary_indices(mb, iface, Val(N))
  left_ghost_axis = _face_axis_index(
    Tuple(left_block.iterators.cell.domain.indices), iface.left, 1
  )
  right_ghost_axis = _face_axis_index(
    Tuple(right_block.iterators.cell.domain.indices), iface.right, 1
  )
  left_ghost = [
    _replace_axis(idx, iface.left.axis, left_ghost_axis, Val(N)) for idx in left_interior
  ]
  right_ghost = [
    _replace_axis(idx, iface.right.axis, right_ghost_axis, Val(N)) for idx in right_interior
  ]

  if length(left_interior) != length(right_interior) ||
    length(left_interior) != length(left_ghost) ||
    length(right_interior) != length(right_ghost)
    throw(
      ArgumentError("Interface cache build failed due to inconsistent face index lengths.")
    )
  end

  A_lr = Vector{SMatrix{N,N,T}}(undef, length(left_interior))
  A_rl = Vector{SMatrix{N,N,T}}(undef, length(left_interior))
  for i in eachindex(left_interior)
    ql = _as_coord_svector(centroid(left_block, left_interior[i]), T, N)
    qr = _as_coord_svector(centroid(right_block, right_interior[i]), T, N)
    A_lr[i] = _basis_transfer_matrix(left_block, ql, right_block, qr, Val(N), T)
    A_rl[i] = _basis_transfer_matrix(right_block, qr, left_block, ql, Val(N), T)
  end

  return InterfaceCache{N,T}(
    left_interior,
    right_interior,
    left_ghost,
    right_ghost,
    A_lr,
    A_rl,
    _grid_state_token(left_block),
    _grid_state_token(right_block),
  )
end

@inline function _replace_axis(
  idx::CartesianIndex{N}, axis::Int, axis_value::Int, ::Val{N}
) where {N}
  vals = idx.I
  return CartesianIndex(ntuple(i -> (i == axis ? axis_value : vals[i]), N))
end

@inline function _basis_transfer_matrix(
  donor::Union{AbstractMappedOrDiscreteGrid,AbstractOrthogonalGrid},
  q_donor::SVector{N,T},
  receiver::Union{AbstractMappedOrDiscreteGrid,AbstractOrthogonalGrid},
  q_receiver::SVector{N,T},
  ::Val{N},
  ::Type{T},
) where {N,T}
  return basis_transfer_matrix(donor, q_donor, receiver, q_receiver)
end
