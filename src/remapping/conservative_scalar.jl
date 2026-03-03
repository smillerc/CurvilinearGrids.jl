@inline function _grid_state_token(grid::AbstractMappedOrDiscreteGrid)
  return UInt(hash(grid.state[]))
end

@inline function _as_coord_svector(x, ::Type{T}, ::Val{N}) where {T,N}
  return SVector{N,T}(ntuple(i -> T(x[i]), N))
end

@inline function _coord_to_cartesian(::CoordinateSystemTrait, q::SVector{N,T}) where {N,T}
  return q
end

@inline function _coord_to_cartesian(::SphericalCS, q::SVector{3,T}) where {T}
  r, theta, phi = q
  st = sin(theta)
  ct = cos(theta)
  sp = sin(phi)
  cp = cos(phi)
  return SVector{3,T}(r * st * cp, r * st * sp, r * ct)
end

@inline function _cartesian_to_coord(::CoordinateSystemTrait, x::SVector{N,T}) where {N,T}
  return x
end
@inline function _cartesian_to_coord(::SphericalCS, x::SVector{3,T}) where {T}
  xx, yy, zz = x
  r = sqrt(xx * xx + yy * yy + zz * zz)
  if r <= eps(T)
    return SVector{3,T}(zero(T), zero(T), zero(T))
  end
  theta = acos(clamp(zz / r, -one(T), one(T)))
  phi = atan(yy, xx)
  if phi < zero(T)
    phi += T(2pi)
  end
  return SVector{3,T}(r, theta, phi)
end

@inline function _cartesian_bounds(
  grid::AbstractMappedOrDiscreteGrid, ::Val{N}, ::Type{T}
) where {N,T}
  mins = MVector{N,T}(ntuple(_ -> typemax(T), N))
  maxs = MVector{N,T}(ntuple(_ -> typemin(T), N))
  cs = coordinate_system(grid)

  for I in grid.iterators.node.domain
    q = SVector{N,T}(ntuple(d -> T(grid.node_coordinates[d][I]), N))
    x = _coord_to_cartesian(cs, q)
    for d in 1:N
      mins[d] = min(mins[d], x[d])
      maxs[d] = max(maxs[d], x[d])
    end
  end

  return SVector{N,T}(mins), SVector{N,T}(maxs)
end

@inline function _in_cartesian_bounds(
  x::SVector{N,T}, xmin::SVector{N,T}, xmax::SVector{N,T}, tol::T
) where {N,T}
  for d in 1:N
    if x[d] < xmin[d] - tol || x[d] > xmax[d] + tol
      return false
    end
  end
  return true
end

@inline function _quadrature_rule(::Type{T}, order::Int) where {T}
  if order == 1
    return T[zero(T)], T[one(T)]
  elseif order == 2
    a = inv(T(2) * sqrt(T(3)))
    return T[-a, a], T[one(T) / 2, one(T) / 2]
  elseif order == 3
    a = sqrt(T(3) / T(5)) / T(2)
    return T[-a, zero(T), a], T[T(5) / 18, T(4) / 9, T(5) / 18]
  end
  throw(
    ArgumentError(
      "Unsupported `quadrature_order=$order`. Supported values are 1, 2, and 3."
    ),
  )
end

@inline function _quadrature_offsets_weights(
  ::Val{1}, points::AbstractVector{T}, weights::AbstractVector{T}
) where {T}
  n = length(points)
  offsets = Vector{SVector{1,T}}(undef, n)
  qweights = Vector{T}(undef, n)
  for i in 1:n
    offsets[i] = SVector{1,T}(points[i])
    qweights[i] = weights[i]
  end
  return offsets, qweights
end

@inline function _quadrature_offsets_weights(
  ::Val{2}, points::AbstractVector{T}, weights::AbstractVector{T}
) where {T}
  n = length(points)
  ns = n * n
  offsets = Vector{SVector{2,T}}(undef, ns)
  qweights = Vector{T}(undef, ns)
  k = 1
  for i in 1:n
    for j in 1:n
      offsets[k] = SVector{2,T}(points[i], points[j])
      qweights[k] = weights[i] * weights[j]
      k += 1
    end
  end
  return offsets, qweights
end

@inline function _quadrature_offsets_weights(
  ::Val{3}, points::AbstractVector{T}, weights::AbstractVector{T}
) where {T}
  n = length(points)
  ns = n * n * n
  offsets = Vector{SVector{3,T}}(undef, ns)
  qweights = Vector{T}(undef, ns)
  k = 1
  for i in 1:n
    for j in 1:n
      for l in 1:n
        offsets[k] = SVector{3,T}(points[i], points[j], points[l])
        qweights[k] = weights[i] * weights[j] * weights[l]
        k += 1
      end
    end
  end
  return offsets, qweights
end

@inline function _source_cell_center(
  source::AbstractMappedOrDiscreteGrid,
  source_global::CartesianIndex{N},
  ::Type{T},
  ::Val{N},
) where {N,T}
  return SVector{N,T}(ntuple(d -> T(source_global.I[d] - source.nhalo) + T(0.5), N))
end

@inline function _destination_local_cell_index(
  destination::AbstractMappedOrDiscreteGrid, xi::SVector{N,T}, tol::T, ::Val{N}
) where {N,T}
  ranges = Tuple(destination.iterators.cell.domain.indices)
  idx = MVector{N,Int}(undef)

  for d in 1:N
    n = length(ranges[d])
    lower_face = T(first(ranges[d]) - destination.nhalo)
    u = xi[d] - lower_face
    if u < -tol || u > T(n) + tol
      return nothing
    end
    i = floor(Int, u) + 1
    idx[d] = clamp(i, 1, n)
  end

  return CartesianIndex(Tuple(idx))
end

@inline function _source_values_vector(
  cache::RemapCache{N,T}, source_field::AbstractArray{<:Real,N}
) where {N,T}
  if size(source_field) != cache.source_shape
    throw(
      ArgumentError(
        "Source field size $(size(source_field)) does not match cache source shape $(cache.source_shape).",
      ),
    )
  end
  return T.(vec(source_field))
end

"""
    validate_remap_cache(cache, source, destination)

Validate that a `RemapCache` still matches the source and destination grids used
to build it.

# Arguments
  - `cache`: Remap cache instance.
  - `source`: Source unified grid.
  - `destination`: Destination unified grid.

# Returns
`true` when validation succeeds; otherwise throws `ArgumentError`.
"""
function validate_remap_cache(
  cache::RemapCache{N},
  source::AbstractMappedOrDiscreteGrid,
  destination::AbstractMappedOrDiscreteGrid,
) where {N}
  if cellsize(source) != cache.source_shape
    throw(
      ArgumentError(
        "Source grid shape $(cellsize(source)) does not match cache source shape $(cache.source_shape).",
      ),
    )
  end
  if cellsize(destination) != cache.destination_shape
    throw(
      ArgumentError(
        "Destination grid shape $(cellsize(destination)) does not match cache destination shape $(cache.destination_shape).",
      ),
    )
  end
  if _grid_state_token(source) != cache.source_state_token
    throw(
      ArgumentError(
        "Source grid state does not match cache state token. Rebuild `RemapCache` after updating the source grid.",
      ),
    )
  end
  if _grid_state_token(destination) != cache.destination_state_token
    throw(
      ArgumentError(
        "Destination grid state does not match cache state token. Rebuild `RemapCache` after updating the destination grid.",
      ),
    )
  end
  return true
end

"""
    build_remap_cache(source, destination; kwargs...)

Build a reusable conservative overlap cache for scalar remapping from `source`
to `destination`.

The cache stores approximate source-to-destination overlap volumes assembled by
sampling each source cell using tensor-product Gauss-Legendre quadrature in the
source computational space and mapping each sample into destination cell space.

# Arguments
  - `source`: Source `MappedGrid` or `DiscreteGrid`.
  - `destination`: Destination `MappedGrid` or `DiscreteGrid`.

# Keywords
  - `quadrature_order`: Per-axis quadrature order (1, 2, or 3). Default: `2`.
  - `inverse_tol`: Inverse-map tolerance for destination coordinate solve. Default: `sqrt(eps(T))`.
  - `inverse_maxiters`: Max Newton iterations for inverse-map solve. Default: `25`.
  - `inverse_bounds_tolerance`: Bounds tolerance for destination inverse-map solve. Default: `sqrt(eps(T))`.
  - `destination_index_tolerance`: Tolerance used when binning destination computational coordinates to destination cells. Default: `sqrt(eps(T))`.
  - `destination_bbox_tolerance`: Cartesian bounding-box prefilter tolerance. Default: `sqrt(eps(T))`.

# Returns
`RemapCache`.
"""
function build_remap_cache(
  source::AbstractMappedOrDiscreteGrid,
  destination::AbstractMappedOrDiscreteGrid;
  quadrature_order::Int=2,
  inverse_tol::Real=sqrt(eps(promote_type(eltype(source), eltype(destination)))),
  inverse_maxiters::Int=25,
  inverse_bounds_tolerance::Real=sqrt(
    eps(promote_type(eltype(source), eltype(destination)))
  ),
  destination_index_tolerance::Real=sqrt(
    eps(promote_type(eltype(source), eltype(destination)))
  ),
  destination_bbox_tolerance::Real=sqrt(
    eps(promote_type(eltype(source), eltype(destination)))
  ),
)
  source_shape = cellsize(source)
  destination_shape = cellsize(destination)
  nsrc_dims = length(source_shape)
  ndst_dims = length(destination_shape)
  if nsrc_dims != ndst_dims
    throw(
      ArgumentError(
        "Source dimension ($nsrc_dims) must match destination dimension ($ndst_dims)."
      ),
    )
  end
  if nsrc_dims < 1 || nsrc_dims > 3
    throw(
      ArgumentError(
        "Conservative scalar remapping currently supports only 1D, 2D, and 3D grids."
      ),
    )
  end

  N = nsrc_dims
  T = promote_type(eltype(source), eltype(destination))
  inv_tol = T(inverse_tol)
  inv_bbox_tol = T(inverse_bounds_tolerance)
  index_tol = T(destination_index_tolerance)
  bbox_tol = T(destination_bbox_tolerance)

  source_volumes = vec(T.(cellvolumes(source)))
  destination_volumes = vec(T.(cellvolumes(destination)))

  points, weights = _quadrature_rule(T, quadrature_order)
  offsets, qweights = _quadrature_offsets_weights(Val(N), points, weights)

  source_linear = LinearIndices(source_shape)
  destination_linear = LinearIndices(destination_shape)
  source_global = source.iterators.global_domain.cell.domain
  source_cs = coordinate_system(source)
  destination_cs = coordinate_system(destination)
  destination_xmin, destination_xmax = _cartesian_bounds(destination, Val(N), T)

  row_idx = Int[]
  col_idx = Int[]
  overlap_values = T[]
  sizehint!(row_idx, length(source_volumes) * max(1, length(qweights) ÷ 2))
  sizehint!(col_idx, length(source_volumes) * max(1, length(qweights) ÷ 2))
  sizehint!(overlap_values, length(source_volumes) * max(1, length(qweights) ÷ 2))

  for source_local in CartesianIndices(source_shape)
    source_lin = source_linear[source_local]
    source_global_idx = source_global[source_local]
    xi_center = _source_cell_center(source, source_global_idx, T, Val(N))
    bins = Dict{Int,T}()

    for q in eachindex(offsets)
      xi_sample = xi_center + offsets[q]
      q_source = _as_coord_svector(coord(source, Tuple(xi_sample)), T, Val(N))
      x_sample = _coord_to_cartesian(source_cs, q_source)
      if !_in_cartesian_bounds(x_sample, destination_xmin, destination_xmax, bbox_tol)
        continue
      end

      q_destination = _cartesian_to_coord(destination_cs, x_sample)
      inverse_result = computational_coordinate(
        destination,
        Tuple(q_destination);
        throw_on_failure=false,
        return_result=true,
        tol=inv_tol,
        maxiters=inverse_maxiters,
        bounds_tolerance=inv_bbox_tol,
      )
      if !inverse_result.converged
        continue
      end

      destination_local = _destination_local_cell_index(
        destination, inverse_result.coordinate, index_tol, Val(N)
      )
      destination_local === nothing && continue
      destination_lin = destination_linear[destination_local]
      bins[destination_lin] = get(bins, destination_lin, zero(T)) + qweights[q]
    end

    source_volume = source_volumes[source_lin]
    for (destination_lin, overlap_fraction) in bins
      overlap_volume = source_volume * overlap_fraction
      if overlap_volume > zero(T)
        push!(row_idx, destination_lin)
        push!(col_idx, source_lin)
        push!(overlap_values, overlap_volume)
      end
    end
  end

  nsrc = length(source_volumes)
  ndst = length(destination_volumes)
  overlap_volume = sparse(row_idx, col_idx, overlap_values, ndst, nsrc)
  source_overlap_volumes = Vector{T}(vec(sum(overlap_volume; dims=1)))
  destination_overlap_volumes = Vector{T}(vec(sum(overlap_volume; dims=2)))

  return RemapCache(
    source_shape,
    destination_shape,
    overlap_volume,
    source_volumes,
    destination_volumes,
    source_overlap_volumes,
    destination_overlap_volumes,
    _grid_state_token(source),
    _grid_state_token(destination),
    quadrature_order,
    length(qweights),
  )
end

"""
    source_overlap_mass(cache, source_field)

Compute the total source mass represented by the overlap region captured in a
`RemapCache`.

# Arguments
  - `cache`: Remap cache instance.
  - `source_field`: Source cell-centered scalar field (cell averages), sized to
    `cache.source_shape`.

# Returns
Overlap mass scalar.
"""
function source_overlap_mass(
  cache::RemapCache{N,T}, source_field::AbstractArray{<:Real,N}
) where {N,T}
  source_values = _source_values_vector(cache, source_field)
  return dot(cache.source_overlap_volumes, source_values)
end

"""
    remap_scalar(cache, source_field; fill_value=0)

Conservatively remap one scalar cell-centered field using a prebuilt `RemapCache`.

# Arguments
  - `cache`: Remap cache instance.
  - `source_field`: Source cell-centered scalar field (cell averages), sized to
    `cache.source_shape`.

# Keywords
  - `fill_value`: Value assigned to destination cells with zero overlap coverage. Default: `0`.

# Returns
Destination remapped scalar field sized to `cache.destination_shape`.
"""
function remap_scalar(
  cache::RemapCache{N,T}, source_field::AbstractArray{<:Real,N}; fill_value::Real=zero(T)
) where {N,T}
  destination_field = Array{T}(undef, cache.destination_shape)
  remap_scalar!(destination_field, cache, source_field; fill_value=fill_value)
  return destination_field
end

"""
    remap_scalar!(destination_field, cache, source_field; fill_value=0)

In-place conservative scalar remap using a prebuilt `RemapCache`.

# Arguments
  - `destination_field`: Destination array sized to `cache.destination_shape`.
  - `cache`: Remap cache instance.
  - `source_field`: Source cell-centered scalar field (cell averages), sized to
    `cache.source_shape`.

# Keywords
  - `fill_value`: Value assigned to destination cells with zero overlap coverage. Default: `0`.

# Returns
`destination_field`.
"""
function remap_scalar!(
  destination_field::AbstractArray{Td,N},
  cache::RemapCache{N,T},
  source_field::AbstractArray{<:Real,N};
  fill_value::Real=zero(T),
) where {N,T,Td<:Real}
  if size(destination_field) != cache.destination_shape
    throw(
      ArgumentError(
        "Destination field size $(size(destination_field)) does not match cache destination shape $(cache.destination_shape).",
      ),
    )
  end

  source_values = _source_values_vector(cache, source_field)
  destination_mass = cache.overlap_volume * source_values
  destination_values = vec(destination_field)
  fill_t = Td(fill_value)

  for i in eachindex(destination_values)
    if cache.destination_overlap_volumes[i] > zero(T) &&
      cache.destination_cell_volumes[i] > zero(T)
      destination_values[i] = Td(destination_mass[i] / cache.destination_cell_volumes[i])
    else
      destination_values[i] = fill_t
    end
  end
  return destination_field
end
