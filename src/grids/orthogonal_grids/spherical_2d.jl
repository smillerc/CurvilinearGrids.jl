function SphericalOrthogonalGrid2D(
  _r::AbstractVector{T},
  _θ::AbstractVector{T},
  nhalo::Int,
  backend;
  halo_coords_included=false,
) where {T<:Real}
  coords, limits, iters, nodedims, celldims = _prepare_nd_coordinates(
    (_r, _θ), nhalo, halo_coords_included
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

  compute_spherical_2d_centroids!(
    centroid_coordinates,
    node_coordinates,
    iters,
    backend,
    halo_coords_included,
  )
  compute_spherical_2d_volumes!(
    cell_volumes,
    node_coordinates,
    iters,
    backend,
    halo_coords_included,
  )
  compute_spherical_2d_face_areas!(
    face_areas,
    node_coordinates,
    iters,
    backend,
    halo_coords_included,
    nhalo,
  )

  return OrthogonalGrid{
    2,
    T,
    SphericalCS,
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

function compute_spherical_2d_centroids!(
  centroids,
  node_coordinates,
  iters,
  backend,
  halo_coords_included,
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_spherical_2d_centroids!(backend)(
    centroids[1],
    centroids[2],
    node_coordinates[1],
    node_coordinates[2],
    domain;
    ndrange=size(domain),
  )
  return nothing
end

function compute_spherical_2d_face_areas!(
  face_areas,
  node_coordinates,
  iters,
  backend,
  halo_coords_included,
  nhalo,
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

  _compute_spherical_2d_radial_face_areas!(backend)(
    face_areas[1],
    node_coordinates[1],
    node_coordinates[2],
    i₊½_domain;
    ndrange=size(i₊½_domain),
  )
  _compute_spherical_2d_theta_face_areas!(backend)(
    face_areas[2],
    node_coordinates[1],
    node_coordinates[2],
    j₊½_domain;
    ndrange=size(j₊½_domain),
  )
  return nothing
end

function compute_spherical_2d_volumes!(
  volumes,
  node_coordinates,
  iters,
  backend,
  halo_coords_included,
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_spherical_2d_volumes!(backend)(
    volumes,
    node_coordinates[1],
    node_coordinates[2],
    domain;
    ndrange=size(domain),
  )
  return nothing
end

@kernel function _compute_spherical_2d_centroids!(rc, θc, rnode, θnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j = I.I

  r₀ = rnode[i]
  r₁ = rnode[i + 1]
  θ₀ = θnode[j]
  θ₁ = θnode[j + 1]

  rc[i] = (3 / 4) * ((r₁^4 - r₀^4) / (r₁^3 - r₀^3))
  θc[j] = acos((cos(θ₀) + cos(θ₁)) / 2)
end

@kernel function _compute_spherical_2d_radial_face_areas!(Aᵢ₊½, r, θ, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j = I.I

  Δμ = cos(θ[j]) - cos(θ[j + 1])
  Aᵢ₊½[I] = 2π * r[i + 1]^2 * Δμ
end

@kernel function _compute_spherical_2d_theta_face_areas!(Aⱼ₊½, r, θ, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j = I.I

  Δr² = r[i + 1]^2 - r[i]^2
  Aⱼ₊½[I] = π * Δr² * sin(θ[j + 1])
end

@kernel function _compute_spherical_2d_volumes!(V, r, θ, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j = I.I

  Δr³ = r[i + 1]^3 - r[i]^3
  Δμ = cos(θ[j]) - cos(θ[j + 1])
  V[I] = (2π / 3) * Δr³ * Δμ
end
