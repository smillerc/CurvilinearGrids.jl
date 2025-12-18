struct SphericalOrthogonalGrid1D{NC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid1D
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

function SphericalOrthogonalGrid1D(
  _r::AbstractVector{T}, nhalo::Int, backend; halo_coords_included=false
) where {T<:Real}
  r, limits, iters, nodedims, celldims = _prepare_1d_coordinates(
    _r, nhalo, halo_coords_included
  )

  node_coordinates = (; r=KernelAbstractions.zeros(backend, T, nodedims...))
  centroid_coordinates = (; r=KernelAbstractions.zeros(backend, T, celldims...))
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (; i₊½=KernelAbstractions.zeros(backend, T, celldims...))

  _populate_1d_nodes!(node_coordinates.r, r, iters, halo_coords_included)

  compute_spherical_1d_centroids!(
    centroid_coordinates, node_coordinates, iters, backend, halo_coords_included
  )
  compute_spherical_1d_volumes!(
    cell_volumes, node_coordinates, iters, backend, halo_coords_included
  )
  compute_spherical_1d_face_areas!(
    face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
  )

  return SphericalOrthogonalGrid1D(
    node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo
  )
end

function compute_spherical_1d_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_spherical_1d_centroids!(backend)(
    centroids.r, node_coordinates.r, domain; ndrange=size(domain)
  )
  return nothing
end

function compute_spherical_1d_face_areas!(
  face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
)
  domain = iters.cell.domain
  if halo_coords_included
    offset = nhalo
    i₊½_domain = expand_upper(expand_lower(domain, 1, offset), 1, offset - 1)
  else
    offset = 1
    i₊½_domain = expand_lower(domain, 1, offset)
  end

  _compute_spherical_1d_face_areas!(backend)(
    face_areas.i₊½, node_coordinates.r, i₊½_domain; ndrange=size(i₊½_domain)
  )
  return nothing
end

function compute_spherical_1d_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_spherical_1d_volumes!(backend)(
    volumes, node_coordinates.r, domain; ndrange=size(domain)
  )
  return nothing
end

@kernel function _compute_spherical_1d_centroids!(rc, rnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I

  r₀ = rnode[i]
  r₁ = rnode[i + 1]

  rc[i] = (3 / 4) * ((r₁^4 - r₀^4) / (r₁^3 - r₀^3))
end

@kernel function _compute_spherical_1d_face_areas!(Aᵢ₊½, r, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, = I.I
  Aᵢ₊½[I] = 4π * r[i + 1]^2
end

@kernel function _compute_spherical_1d_volumes!(V, r, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I
  V[I] = (4 / 3) * π * (r[i + 1]^3 - r[i]^3)
end