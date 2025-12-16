struct AxisymmetricOrthogonalGrid2D{NC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid2D
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

function AxisymmetricOrthogonalGrid2D(
  _r::AbstractVector{T},
  _z::AbstractVector{T},
  nhalo::Int,
  backend;
  halo_coords_included=false,
) where {T<:Real}
  coords, limits, iters, nodedims, celldims = _prepare_nd_coordinates(
    (_r, _z), nhalo, halo_coords_included
  )

  node_coordinates = (;
    r=KernelAbstractions.zeros(backend, T, nodedims[1]),
    z=KernelAbstractions.zeros(backend, T, nodedims[2]),
  )
  centroid_coordinates = (;
    r=KernelAbstractions.zeros(backend, T, celldims[1]),
    z=KernelAbstractions.zeros(backend, T, celldims[2]),
  )
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (;
    i₊½=KernelAbstractions.zeros(backend, T, celldims...),
    j₊½=KernelAbstractions.zeros(backend, T, celldims...),
  )

  _populate_2d_nodes!(node_coordinates, coords, iters, halo_coords_included)

  compute_axisymmetric_centroids!(
    centroid_coordinates, node_coordinates, iters, backend, halo_coords_included
  )
  compute_axisymmetric_volumes!(
    cell_volumes, node_coordinates, iters, backend, halo_coords_included
  )
  compute_axisymmetric_face_areas!(
    face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
  )

  return AxisymmetricOrthogonalGrid2D(
    node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo
  )
end

function _populate_2d_nodes!(storage, coords, iters, halo_coords_included)
  r, z = coords
  if halo_coords_included
    copy!(storage.r, r)
    copy!(storage.z, z)
  else
    @views begin
      copy!(storage.r[iters.node.domain.indices[1]], r)
      copy!(storage.z[iters.node.domain.indices[2]], z)
    end
  end
  return nothing
end

function compute_axisymmetric_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_axisymmetric_centroids!(backend)(
    centroids.r,
    centroids.z,
    node_coordinates.r,
    node_coordinates.z,
    domain;
    ndrange=size(domain),
  )
  return nothing
end

function compute_axisymmetric_face_areas!(
  face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
)
  domain = iters.cell.domain
  if halo_coords_included
    offset = +nhalo
    i₊½_domain = expand_upper(expand_lower(domain, 1, offset), 1, offset - 1)
    j₊½_domain = expand_upper(expand_lower(domain, 2, offset), 2, offset - 1)
  else
    offset = +1
    i₊½_domain = expand_lower(domain, 1, offset)
    j₊½_domain = expand_lower(domain, 2, offset)
  end

  _compute_axisymmetric_radial_face_areas!(backend)(
    face_areas.i₊½,
    node_coordinates.r,
    node_coordinates.z,
    i₊½_domain;
    ndrange=size(i₊½_domain),
  )

  _compute_axisymmetric_axial_face_areas!(backend)(
    face_areas.j₊½, node_coordinates.r, j₊½_domain; ndrange=size(j₊½_domain)
  )
  return nothing
end

function compute_axisymmetric_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_axisymmetric_volumes!(backend)(
    volumes, node_coordinates.r, node_coordinates.z, domain; ndrange=size(domain)
  )
  return nothing
end

@kernel function _compute_axisymmetric_centroids!(rc, zc, rnode, znode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j = I.I

  r₀ = rnode[i]
  r₁ = rnode[i + 1]
  z₀ = znode[j]
  z₁ = znode[j + 1]

  rc[i] = (2 / 3) * ((r₁^3 - r₀^3) / (r₁^2 - r₀^2))
  zc[j] = (z₀ + z₁) / 2
end

@kernel function _compute_axisymmetric_radial_face_areas!(Aᵢ₊½, r, z, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j = I.I

  Δz = z[j + 1] - z[j]
  Aᵢ₊½[I] = 2π * r[i + 1] * Δz
end

@kernel function _compute_axisymmetric_axial_face_areas!(Aⱼ₊½, r, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j = I.I

  Aⱼ₊½[I] = π * (r[i + 1]^2 - r[i]^2)
end

@kernel function _compute_axisymmetric_volumes!(V, r, z, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j = I.I

  Δr² = r[i + 1]^2 - r[i]^2
  Δz = z[j + 1] - z[j]

  V[I] = π * Δr² * Δz
end