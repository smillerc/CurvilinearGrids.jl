function CartesianOrthogonalGrid3D(
  _x::AbstractVector{T},
  _y::AbstractVector{T},
  _z::AbstractVector{T},
  nhalo::Int,
  backend;
  halo_coords_included=false,
) where {T<:Real}
  coords, limits, iters, nodedims, celldims = _prepare_nd_coordinates(
    (_x, _y, _z), nhalo, halo_coords_included
  )

  node_coordinates = (
    KernelAbstractions.zeros(backend, T, nodedims[1]),
    KernelAbstractions.zeros(backend, T, nodedims[2]),
    KernelAbstractions.zeros(backend, T, nodedims[3]),
  )
  centroid_coordinates = (
    KernelAbstractions.zeros(backend, T, celldims[1]),
    KernelAbstractions.zeros(backend, T, celldims[2]),
    KernelAbstractions.zeros(backend, T, celldims[3]),
  )
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (
    KernelAbstractions.zeros(backend, T, celldims...),
    KernelAbstractions.zeros(backend, T, celldims...),
    KernelAbstractions.zeros(backend, T, celldims...),
  )

  _populate_3d_nodes!(node_coordinates, coords, iters)

  compute_cartesian_3d_centroids!(
    centroid_coordinates, node_coordinates, iters, backend, true
  )
  compute_cartesian_3d_volumes!(
    cell_volumes, node_coordinates, iters, backend, true
  )
  compute_cartesian_3d_face_areas!(
    face_areas, node_coordinates, iters, backend, true, nhalo
  )

  return OrthogonalGrid{
    3,
    T,
    CartesianCS,
    typeof(node_coordinates),
    typeof(centroid_coordinates),
    typeof(cell_volumes),
    typeof(iters),
    typeof(limits),
    typeof(face_areas),
  }(node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo)
end

function _populate_3d_nodes!(storage, coords, iters)
  x, y, z = coords
  if length(x) == length(storage[1]) &&
     length(y) == length(storage[2]) &&
     length(z) == length(storage[3])
    copy!(storage[1], x)
    copy!(storage[2], y)
    copy!(storage[3], z)
  else
    @views begin
      copy!(storage[1][iters.node.domain.indices[1]], x)
      copy!(storage[2][iters.node.domain.indices[2]], y)
      copy!(storage[3][iters.node.domain.indices[3]], z)
    end
  end
  return nothing
end

function compute_cartesian_3d_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain
  _compute_cartesian_3d_centroids!(backend)(
    centroids[1],
    centroids[2],
    centroids[3],
    node_coordinates[1],
    node_coordinates[2],
    node_coordinates[3],
    domain;
    ndrange=size(domain),
  )
  return nothing
end

function compute_cartesian_3d_face_areas!(
  face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
)
  domain = iters.cell.domain
  if nhalo == 0
    i₊½_domain = domain
    j₊½_domain = domain
    k₊½_domain = domain
  elseif halo_coords_included
    offset = nhalo
    i₊½_domain = expand_upper(expand_lower(domain, 1, offset), 1, offset - 1)
    j₊½_domain = expand_upper(expand_lower(domain, 2, offset), 2, offset - 1)
    k₊½_domain = expand_upper(expand_lower(domain, 3, offset), 3, offset - 1)
  else
    offset = 1
    i₊½_domain = expand_lower(domain, 1, offset)
    j₊½_domain = expand_lower(domain, 2, offset)
    k₊½_domain = expand_lower(domain, 3, offset)
  end

  _compute_cartesian_3d_face_areas_i!(backend)(
    face_areas[1], node_coordinates[2], node_coordinates[3], i₊½_domain; ndrange=size(i₊½_domain)
  )
  _compute_cartesian_3d_face_areas_j!(backend)(
    face_areas[2], node_coordinates[1], node_coordinates[3], j₊½_domain; ndrange=size(j₊½_domain)
  )
  _compute_cartesian_3d_face_areas_k!(backend)(
    face_areas[3], node_coordinates[1], node_coordinates[2], k₊½_domain; ndrange=size(k₊½_domain)
  )
  return nothing
end

function compute_cartesian_3d_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain
  _compute_cartesian_3d_volumes!(backend)(
    volumes,
    node_coordinates[1],
    node_coordinates[2],
    node_coordinates[3],
    domain;
    ndrange=size(domain),
  )
  return nothing
end

@kernel function _compute_cartesian_3d_centroids!(xc, yc, zc, xnode, ynode, znode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j, k = I.I
  xc[i] = (xnode[i + 1] + xnode[i]) / 2
  yc[j] = (ynode[j + 1] + ynode[j]) / 2
  zc[k] = (znode[k + 1] + znode[k]) / 2
end

@kernel function _compute_cartesian_3d_face_areas_i!(Aᵢ₊½, y, z, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I
  _ = i
  Aᵢ₊½[I] = (y[j + 1] - y[j]) * (z[k + 1] - z[k])
end

@kernel function _compute_cartesian_3d_face_areas_j!(Aⱼ₊½, x, z, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I
  _ = j
  Aⱼ₊½[I] = (x[i + 1] - x[i]) * (z[k + 1] - z[k])
end

@kernel function _compute_cartesian_3d_face_areas_k!(Aₖ₊½, x, y, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I
  _ = k
  Aₖ₊½[I] = (x[i + 1] - x[i]) * (y[j + 1] - y[j])
end

@kernel function _compute_cartesian_3d_volumes!(V, x, y, z, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j, k = I.I
  V[I] = (x[i + 1] - x[i]) * (y[j + 1] - y[j]) * (z[k + 1] - z[k])
end
