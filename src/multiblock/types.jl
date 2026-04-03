"""
    BlockFace{N}

Identify one face of one `N`-dimensional block.

# Fields
  - `block_id`: 1-based block index into `MultiBlockMesh.blocks`.
  - `axis`: Normal axis index in `1:N`.
  - `side`: Face side symbol (`:min` or `:max`).
"""
struct BlockFace{N}
  block_id::Int
  axis::Int
  side::Symbol

  function BlockFace{N}(block_id::Int, axis::Int, side::Symbol) where {N}
    if axis < 1 || axis > N
      throw(ArgumentError("Invalid `axis=$axis` for `N=$N`."))
    end
    if side !== :min && side !== :max
      throw(ArgumentError("Invalid `side=$side`. Expected `:min` or `:max`."))
    end
    return new{N}(block_id, axis, side)
  end
end

"""
    BlockInterface{N,T}

One 1-to-1 interface between two block faces.

# Fields
  - `left`: Left face descriptor.
  - `right`: Right face descriptor.
  - `permutation`: Tangential-axis permutation mapping left -> right.
  - `flips`: Tangential-axis reversal flags (after permutation).
  - `tolerance`: Abutting geometry tolerance.
"""
struct BlockInterface{N,T,P<:Tuple,F<:Tuple}
  left::BlockFace{N}
  right::BlockFace{N}
  permutation::P
  flips::F
  tolerance::T
end

function BlockInterface(
  left::BlockFace{N},
  right::BlockFace{N},
  permutation::Tuple,
  flips::Tuple;
  tolerance::T=T(1e-10),
) where {N,T<:Real}
  if left == right
    throw(ArgumentError("Interface cannot connect a face to itself."))
  end
  if length(permutation) != N - 1
    throw(
      ArgumentError(
        "Invalid `permutation` length $(length(permutation)); expected $(N - 1)."
      ),
    )
  end
  if length(flips) != N - 1
    throw(ArgumentError("Invalid `flips` length $(length(flips)); expected $(N - 1)."))
  end
  if !all(x -> x isa Int, permutation)
    throw(ArgumentError("`permutation` entries must be integers."))
  end
  if !all(x -> x isa Bool, flips)
    throw(ArgumentError("`flips` entries must be booleans."))
  end
  if sort(collect(permutation)) != collect(1:(N - 1))
    throw(
      ArgumentError(
        "Invalid tangential `permutation`. Expected a permutation of `1:$(N - 1)`."
      ),
    )
  end
  if tolerance < 0
    throw(ArgumentError("`tolerance` must be non-negative."))
  end
  perm = Tuple(Int(x) for x in permutation)
  flip = Tuple(Bool(x) for x in flips)
  return BlockInterface{N,T,typeof(perm),typeof(flip)}(left, right, perm, flip, tolerance)
end

"""
    InterfaceCache{N,T}

Cached point-pair connectivity and transforms for one interface.

# Fields
  - `left_interior`: Left interior boundary-cell indices.
  - `right_interior`: Right interior boundary-cell indices.
  - `left_ghost`: Left ghost-cell indices adjacent to interface.
  - `right_ghost`: Right ghost-cell indices adjacent to interface.
  - `left_to_right_transform`: Vector component transforms left -> right.
  - `right_to_left_transform`: Vector component transforms right -> left.
  - `left_state_token`: Left block state token at cache build time.
  - `right_state_token`: Right block state token at cache build time.
"""
struct InterfaceCache{N,T}
  left_interior::Vector{CartesianIndex{N}}
  right_interior::Vector{CartesianIndex{N}}
  left_ghost::Vector{CartesianIndex{N}}
  right_ghost::Vector{CartesianIndex{N}}
  left_to_right_transform::Vector{SMatrix{N,N,T}}
  right_to_left_transform::Vector{SMatrix{N,N,T}}
  left_state_token::UInt
  right_state_token::UInt
end

"""
    MultiBlockCache{C}

Container for per-interface cache payloads.

# Fields
  - `interfaces`: One optional cache per interface.
"""
mutable struct MultiBlockCache{C}
  interfaces::Vector{Union{Nothing,C}}
end

"""
    MultiBlockMesh{N,T,B,I,C}

Topological and geometric connectivity for multiple unified blocks.

# Fields
  - `blocks`: Tuple of unified blocks (`MappedGrid` or `DiscreteGrid`).
  - `interfaces`: Vector of 1-to-1 interface declarations.
  - `cache`: Interface mapping cache container.
"""
struct MultiBlockMesh{N,T,B,I,C}
  blocks::B
  interfaces::I
  cache::C
end

@inline _default_permutation(::Val{N}) where {N} = ntuple(i -> i, N - 1)
@inline _default_flips(::Val{N}) where {N} = ntuple(_ -> false, N - 1)

@inline _normalize_face_symbol(face::Symbol) = Symbol(lowercase(String(face)))

function _face_symbol_to_axis_side(::Val{N}, face::Symbol) where {N}
  f = _normalize_face_symbol(face)
  axis, side = if f in (:ilo, :imin, :xlo, :xmin)
    (1, :min)
  elseif f in (:ihi, :imax, :xhi, :xmax)
    (1, :max)
  elseif f in (:jlo, :jmin, :ylo, :ymin)
    (2, :min)
  elseif f in (:jhi, :jmax, :yhi, :ymax)
    (2, :max)
  elseif f in (:klo, :kmin, :zlo, :zmin)
    (3, :min)
  elseif f in (:khi, :kmax, :zhi, :zmax)
    (3, :max)
  else
    throw(
      ArgumentError(
        "Invalid face symbol `$face`. Expected one of `:ilo/:ihi`, `:jlo/:jhi`, or `:klo/:khi` for higher dimensions.",
      ),
    )
  end

  if axis > N
    throw(
      ArgumentError(
        "Face symbol `$face` maps to axis $axis, which is invalid for `N=$N` dimensions."
      ),
    )
  end
  return axis, side
end

function _find_block_id(blocks::Tuple, mesh)
  for i in eachindex(blocks)
    blocks[i] === mesh && return i
  end
  for i in eachindex(blocks)
    isequal(blocks[i], mesh) && return i
  end
  throw(ArgumentError("Interface references a mesh not present in `blocks`."))
end

function _parse_interface_endpoint(
  endpoint, blocks::Tuple, ::Val{N}, iface_id, side_name
) where {N}
  if !(endpoint isa Tuple) || length(endpoint) != 2
    throw(
      ArgumentError(
        "Interface[$iface_id] $side_name endpoint must be `(mesh, :face_symbol)`, got `$endpoint`.",
      ),
    )
  end

  mesh, face_symbol = endpoint
  if !(face_symbol isa Symbol)
    throw(
      ArgumentError(
        "Interface[$iface_id] $side_name endpoint face must be a `Symbol`, got `$(typeof(face_symbol))`.",
      ),
    )
  end

  block_id = _find_block_id(blocks, mesh)
  axis, side = _face_symbol_to_axis_side(Val(N), face_symbol)
  return BlockFace{N}(block_id, axis, side)
end

function _parse_interface_spec(
  spec, blocks::Tuple, ::Val{N}, ::Type{T}, iface_id::Int, default_tolerance::T
) where {N,T}
  if spec isa BlockInterface{N}
    return BlockInterface(
      spec.left, spec.right, spec.permutation, spec.flips; tolerance=T(spec.tolerance)
    )
  end
  if !(spec isa Pair)
    throw(
      ArgumentError(
        "Interface[$iface_id] must be either `BlockInterface{$N}` or `(mesh, :face) => (mesh, :face)`.",
      ),
    )
  end

  left = _parse_interface_endpoint(spec.first, blocks, Val(N), iface_id, :left)
  right = _parse_interface_endpoint(spec.second, blocks, Val(N), iface_id, :right)
  return BlockInterface(
    left,
    right,
    _default_permutation(Val(N)),
    _default_flips(Val(N));
    tolerance=default_tolerance,
  )
end

function _collect_interface_specs(interfaces)
  if !(interfaces isa Tuple)
    throw(
      ArgumentError(
        "`interfaces` must be a tuple, e.g. `((mesh_1, :ihi) => (mesh_2, :ilo),)`."
      ),
    )
  end
  return collect(interfaces)
end

"""
    MultiBlockMesh(blocks, interfaces; validate=true, build_cache=true, tolerance=1e-10)

Construct a multi-block mesh and optionally validate/build interface caches.

# Arguments
  - `blocks`: Collection of `MappedGrid`/`DiscreteGrid`/`OrthogonalGrid` blocks.
  - `interfaces`: Tuple of interface declarations, either:
    - `BlockInterface{N}` entries, or
    - mesh references using `(mesh, :face_symbol) => (mesh, :face_symbol)`.

# Keywords
  - `validate`: Run topology/abutting validation. Default: `true`.
  - `build_cache`: Build interface caches. Default: `true`.
  - `tolerance`: Default abutting tolerance for pair-style interface specs.

# Returns
`MultiBlockMesh`.
"""
function MultiBlockMesh(
  blocks, interfaces; validate::Bool=true, build_cache::Bool=true, tolerance::Real=1e-10
)
  block_tuple = Tuple(blocks)
  isempty(block_tuple) && throw(ArgumentError("`blocks` cannot be empty."))

  first_block = first(block_tuple)
  if !(
    first_block isa AbstractMappedOrDiscreteGrid || first_block isa AbstractOrthogonalGrid
  )
    throw(
      ArgumentError(
        "`blocks` must contain only `MappedGrid`, `DiscreteGrid`, or `OrthogonalGrid`."
      ),
    )
  end
  N = ndims(first_block.iterators.cell.full)
  T = eltype(first_block)

  for (i, block) in enumerate(block_tuple)
    if !(block isa AbstractMappedOrDiscreteGrid || block isa AbstractOrthogonalGrid)
      throw(
        ArgumentError(
          "Block $i has unsupported type $(typeof(block)); expected `MappedGrid`, `DiscreteGrid`, or `OrthogonalGrid`.",
        ),
      )
    end
    nb = ndims(block.iterators.cell.full)
    if nb != N
      throw(ArgumentError("Block $i has dimension $nb but expected $N."))
    end
  end

  iface_specs = _collect_interface_specs(interfaces)
  tol = T(tolerance)
  iface_vec = if isempty(iface_specs)
    BlockInterface{N,T,NTuple{N - 1,Int},NTuple{N - 1,Bool}}[]
  else
    [
      _parse_interface_spec(spec, block_tuple, Val(N), T, iface_id, tol) for
      (iface_id, spec) in pairs(iface_specs)
    ]
  end

  cache = MultiBlockCache{InterfaceCache{N,T}}(
    Union{Nothing,InterfaceCache{N,T}}[nothing for _ in eachindex(iface_vec)]
  )

  mb = MultiBlockMesh{N,T,typeof(block_tuple),typeof(iface_vec),typeof(cache)}(
    block_tuple, iface_vec, cache
  )

  if validate
    validate_multiblock!(mb)
  end
  if build_cache
    build_interface_caches!(mb)
  end
  return mb
end

@inline _face_key(face::BlockFace{N}) where {N} = (face.block_id, face.axis, face.side)

@inline function _grid_state_token(grid::AbstractMappedOrDiscreteGrid)
  return UInt(hash(grid.state[]))
end

@inline function _grid_state_token(grid::AbstractOrthogonalGrid)
  return UInt(hash(grid.node_coordinates))
end
