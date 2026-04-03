@inline _axisymmetric_rotational_axis(::Symbol) = throw(
  ArgumentError("rotational_axis must be `:x` or `:y`")
)
@inline _axisymmetric_rotational_axis(::Val{:x}) = :x
@inline _axisymmetric_rotational_axis(::Val{:y}) = :y

function AxisymmetricOrthogonalGrid2D(
  _coord1::AbstractVector{T},
  _coord2::AbstractVector{T},
  nhalo::Int,
  backend;
  halo_coords_included=false,
  rotational_axis::Symbol=:y,
) where {T<:Real}
  rotational_axis = _axisymmetric_rotational_axis(Val(rotational_axis))
  coords, limits, iters, nodedims, celldims = _prepare_nd_coordinates(
    (_coord1, _coord2), nhalo, halo_coords_included
  )

  node_coordinates = (
    KernelAbstractions.zeros(backend, T, nodedims[1]),
    KernelAbstractions.zeros(backend, T, nodedims[2]),
  )
  centroid_coordinates = (
    KernelAbstractions.zeros(backend, T, celldims[1]),
    KernelAbstractions.zeros(backend, T, celldims[2]),
  )
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (
    KernelAbstractions.zeros(backend, T, celldims...),
    KernelAbstractions.zeros(backend, T, celldims...),
  )

  _populate_2d_nodes!(node_coordinates, coords, iters)

  compute_axisymmetric_centroids!(
    centroid_coordinates,
    node_coordinates,
    iters,
    backend,
    halo_coords_included,
    rotational_axis,
  )
  compute_axisymmetric_volumes!(
    cell_volumes, node_coordinates, iters, backend, halo_coords_included, rotational_axis
  )
  compute_axisymmetric_face_areas!(
    face_areas,
    node_coordinates,
    iters,
    backend,
    halo_coords_included,
    nhalo,
    rotational_axis,
  )

  coordinate_system = rotational_axis === :x ? AxisymmetricCS{:x}() : AxisymmetricCS{:y}()
  return OrthogonalGrid{
    2,
    T,
    typeof(coordinate_system),
    typeof(node_coordinates),
    typeof(centroid_coordinates),
    typeof(cell_volumes),
    typeof(iters),
    typeof(limits),
    typeof(face_areas),
  }(
    node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo
  )
end

function _populate_2d_nodes!(storage, coords, iters)
  r, z = coords
  if length(r) == length(storage[1]) && length(z) == length(storage[2])
    copy!(storage[1], r)
    copy!(storage[2], z)
  else
    @views begin
      copy!(storage[1][iters.node.domain.indices[1]], r)
      copy!(storage[2][iters.node.domain.indices[2]], z)
    end
  end
  return nothing
end

function compute_axisymmetric_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included, rotational_axis
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain
  radial_dim = rotational_axis === :x ? 2 : 1

  _compute_axisymmetric_centroids!(backend)(
    centroids[1],
    centroids[2],
    node_coordinates[1],
    node_coordinates[2],
    domain,
    Val(radial_dim);
    ndrange=size(domain),
  )
  return nothing
end

function compute_axisymmetric_face_areas!(
  face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo, rotational_axis
)
  domain = iters.cell.domain
  if nhalo == 0
    i₊½_domain = domain
    j₊½_domain = domain
  elseif halo_coords_included
    offset = +nhalo
    i₊½_domain = expand_upper(expand_lower(domain, 1, offset), 1, offset - 1)
    j₊½_domain = expand_upper(expand_lower(domain, 2, offset), 2, offset - 1)
  else
    offset = +1
    i₊½_domain = expand_lower(domain, 1, offset)
    j₊½_domain = expand_lower(domain, 2, offset)
  end
  radial_dim = rotational_axis === :x ? 2 : 1

  _compute_axisymmetric_radial_face_areas!(backend)(
    face_areas[1],
    node_coordinates[1],
    node_coordinates[2],
    i₊½_domain,
    Val(radial_dim);
    ndrange=size(i₊½_domain),
  )

  _compute_axisymmetric_axial_face_areas!(backend)(
    face_areas[2],
    node_coordinates[1],
    node_coordinates[2],
    j₊½_domain,
    Val(radial_dim);
    ndrange=size(j₊½_domain),
  )
  return nothing
end

function compute_axisymmetric_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included, rotational_axis
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain
  radial_dim = rotational_axis === :x ? 2 : 1

  _compute_axisymmetric_volumes!(backend)(
    volumes,
    node_coordinates[1],
    node_coordinates[2],
    domain,
    Val(radial_dim);
    ndrange=size(domain),
  )
  return nothing
end

@inline _axisymmetric_midpoint(a, b) = (a + b) / 2
@inline _axisymmetric_radial_centroid(a, b) = (2 / 3) * ((b^3 - a^3) / (b^2 - a^2))

@kernel function _compute_axisymmetric_centroids!(
  c1, c2, n1, n2, cell_domain, radial_dim::Val{RadialDim}
) where {RadialDim}
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j = I.I

  if RadialDim == 1
    c1[i] = _axisymmetric_radial_centroid(n1[i], n1[i + 1])
    c2[j] = _axisymmetric_midpoint(n2[j], n2[j + 1])
  else
    c1[i] = _axisymmetric_midpoint(n1[i], n1[i + 1])
    c2[j] = _axisymmetric_radial_centroid(n2[j], n2[j + 1])
  end
end

@kernel function _compute_axisymmetric_radial_face_areas!(
  Aᵢ₊½, n1, n2, domain, radial_dim::Val{RadialDim}
) where {RadialDim}
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j = I.I

  if RadialDim == 1
    Aᵢ₊½[I] = 2π * n1[i + 1] * (n2[j + 1] - n2[j])
  else
    Aᵢ₊½[I] = π * (n2[j + 1]^2 - n2[j]^2)
  end
end

@kernel function _compute_axisymmetric_axial_face_areas!(
  Aⱼ₊½, n1, n2, domain, radial_dim::Val{RadialDim}
) where {RadialDim}
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j = I.I

  if RadialDim == 1
    Aⱼ₊½[I] = π * (n1[i + 1]^2 - n1[i]^2)
  else
    Aⱼ₊½[I] = 2π * n2[j + 1] * (n1[i + 1] - n1[i])
  end
end

@kernel function _compute_axisymmetric_volumes!(
  V, n1, n2, cell_domain, radial_dim::Val{RadialDim}
) where {RadialDim}
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j = I.I

  if RadialDim == 1
    V[I] = π * (n1[i + 1]^2 - n1[i]^2) * (n2[j + 1] - n2[j])
  else
    V[I] = π * (n2[j + 1]^2 - n2[j]^2) * (n1[i + 1] - n1[i])
  end
end
