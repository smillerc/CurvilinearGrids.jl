struct CartesianOrthogonalGrid1D{NC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid1D
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

struct CylindricalOrthogonalGrid1D{NC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid1D
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

struct SphericalOrthogonalGrid1D{NC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid1D
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

struct AxisymmetricOrthogonalGrid2D{NC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid2D
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

function CartesianOrthogonalGrid1D(
  _x::AbstractVector{T}, nhalo::Int, backend; halo_coords_included=false
) where {T<:Real}
  x, limits, iters, nodedims, celldims =
    _prepare_1d_coordinates(_x, nhalo, halo_coords_included)

  node_coordinates = (; x=KernelAbstractions.zeros(backend, T, nodedims...))
  centroid_coordinates = (; x=KernelAbstractions.zeros(backend, T, celldims...))
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (; i₊½=KernelAbstractions.zeros(backend, T, celldims...))

  _populate_1d_nodes!(node_coordinates.x, x, iters, halo_coords_included)

  compute_cartesian_centroids!(centroid_coordinates, node_coordinates, iters, backend, halo_coords_included)
  compute_cartesian_volumes!(cell_volumes, node_coordinates, iters, backend, halo_coords_included)
  compute_cartesian_face_areas!(
    face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
  )

  return CartesianOrthogonalGrid1D(
    node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo
  )
end

function CylindricalOrthogonalGrid1D(
  _r::AbstractVector{T}, nhalo::Int, backend; halo_coords_included=false
) where {T<:Real}
  r, limits, iters, nodedims, celldims =
    _prepare_1d_coordinates(_r, nhalo, halo_coords_included)

  node_coordinates = (; r=KernelAbstractions.zeros(backend, T, nodedims...))
  centroid_coordinates = (; r=KernelAbstractions.zeros(backend, T, celldims...))
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (; i₊½=KernelAbstractions.zeros(backend, T, celldims...))

  _populate_1d_nodes!(node_coordinates.r, r, iters, halo_coords_included)

  compute_cylindrical_centroids!(
    centroid_coordinates, node_coordinates, iters, backend, halo_coords_included
  )
  compute_cylindrical_volumes!(
    cell_volumes, node_coordinates, iters, backend, halo_coords_included
  )
  compute_cylindrical_face_areas!(
    face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
  )

  return CylindricalOrthogonalGrid1D(
    node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo
  )
end

function SphericalOrthogonalGrid1D(
  _r::AbstractVector{T}, nhalo::Int, backend; halo_coords_included=false
) where {T<:Real}
  r, limits, iters, nodedims, celldims =
    _prepare_1d_coordinates(_r, nhalo, halo_coords_included)

  node_coordinates = (; r=KernelAbstractions.zeros(backend, T, nodedims...))
  centroid_coordinates = (; r=KernelAbstractions.zeros(backend, T, celldims...))
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (; i₊½=KernelAbstractions.zeros(backend, T, celldims...))

  _populate_1d_nodes!(node_coordinates.r, r, iters, halo_coords_included)

  compute_spherical_centroids!(
    centroid_coordinates, node_coordinates, iters, backend, halo_coords_included
  )
  compute_spherical_volumes!(
    cell_volumes, node_coordinates, iters, backend, halo_coords_included
  )
  compute_spherical_face_areas!(
    face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
  )

  return SphericalOrthogonalGrid1D(
    node_coordinates, centroid_coordinates, cell_volumes, iters, limits, face_areas, nhalo
  )
end

function AxisymmetricOrthogonalGrid2D(
  _r::AbstractVector{T},
  _z::AbstractVector{T},
  nhalo::Int,
  backend;
  halo_coords_included=false,
) where {T<:Real}
  coords, limits, iters, nodedims, celldims =
    _prepare_nd_coordinates((_r, _z), nhalo, halo_coords_included)

  node_coordinates = (
    ; r=KernelAbstractions.zeros(backend, T, nodedims[1]), z=KernelAbstractions.zeros(backend, T, nodedims[2])
  )
  centroid_coordinates = (
    ; r=KernelAbstractions.zeros(backend, T, celldims[1]), z=KernelAbstractions.zeros(backend, T, celldims[2])
  )
  cell_volumes = KernelAbstractions.zeros(backend, T, celldims...)
  face_areas = (
    ; i₊½=KernelAbstractions.zeros(backend, T, celldims...), j₊½=KernelAbstractions.zeros(backend, T, celldims...)
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

function _prepare_1d_coordinates(_x, nhalo, halo_coords_included)
  if !halo_coords_included
    x = pad_with_halo(_x, nhalo)
    halo_coords_included = true
  else
    x = _x
  end

  limits, iters = get_iterators(size(x), halo_coords_included, nhalo)
  celldims = size(iters.cell.full)
  nodedims = size(iters.node.full)

  return x, limits, iters, nodedims, celldims
end

function _prepare_nd_coordinates(_coords::NTuple{2}, nhalo, halo_coords_included)
  if !halo_coords_included
    coords = map(c -> pad_with_halo(c, nhalo), _coords)
    halo_coords_included = true
  else
    coords = _coords
  end

  nodedims = map(length, coords) |> Tuple
  limits, iters = get_iterators(nodedims, halo_coords_included, nhalo)
  celldims = size(iters.cell.full)

  return coords, limits, iters, nodedims, celldims
end

function _populate_1d_nodes!(storage, coords, iters, halo_coords_included)
  if halo_coords_included
    copy!(storage, coords)
  else
    @views copy!(storage[iters.node.domain], coords)
  end
  return nothing
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

function compute_cartesian_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_cartesian_centroids!(backend)(
    centroids.x, node_coordinates.x, domain; ndrange=size(domain)
  )
  return nothing
end

function compute_cartesian_face_areas!(
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

  _compute_cartesian_face_areas!(backend)(
    face_areas.i₊½, node_coordinates.x, i₊½_domain; ndrange=size(i₊½_domain)
  )
  return nothing
end

function compute_cartesian_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_cartesian_volumes!(backend)(
    volumes, node_coordinates.x, domain; ndrange=size(domain)
  )
  return nothing
end

function compute_cylindrical_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_cylindrical_centroids!(backend)(
    centroids.r, node_coordinates.r, domain; ndrange=size(domain)
  )
  return nothing
end

function compute_cylindrical_face_areas!(
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

  _compute_cylindrical_face_areas!(backend)(
    face_areas.i₊½, node_coordinates.r, i₊½_domain; ndrange=size(i₊½_domain)
  )
  return nothing
end

function compute_cylindrical_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_cylindrical_volumes!(backend)(
    volumes, node_coordinates.r, domain; ndrange=size(domain)
  )
  return nothing
end

function compute_spherical_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_spherical_centroids!(backend)(
    centroids.r, node_coordinates.r, domain; ndrange=size(domain)
  )
  return nothing
end

function compute_spherical_face_areas!(
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

  _compute_spherical_face_areas!(backend)(
    face_areas.i₊½, node_coordinates.r, i₊½_domain; ndrange=size(i₊½_domain)
  )
  return nothing
end

function compute_spherical_volumes!(
  volumes, node_coordinates, iters, backend, halo_coords_included
)
  domain = halo_coords_included ? iters.cell.full : iters.cell.domain

  _compute_spherical_volumes!(backend)(
    volumes, node_coordinates.r, domain; ndrange=size(domain)
  )
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
    face_areas.j₊½,
    node_coordinates.r,
    j₊½_domain;
    ndrange=size(j₊½_domain),
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

@kernel function _compute_cartesian_centroids!(xc, xnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I
  xc[i] = (xnode[i + 1] + xnode[i]) / 2
end

@kernel function _compute_cartesian_face_areas!(Aᵢ₊½, x, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, = I.I
  Aᵢ₊½[I] = one(x[1])
end

@kernel function _compute_cartesian_volumes!(V, x, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I
  V[I] = x[i + 1] - x[i]
end

@kernel function _compute_cylindrical_centroids!(rc, rnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I

  r₀ = rnode[i]
  r₁ = rnode[i + 1]

  rc[i] = (2 / 3) * ((r₁^3 - r₀^3) / (r₁^2 - r₀^2))
end

@kernel function _compute_cylindrical_face_areas!(Aᵢ₊½, r, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, = I.I
  Aᵢ₊½[I] = 2π * r[i + 1]
end

@kernel function _compute_cylindrical_volumes!(V, r, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I
  V[I] = π * (r[i + 1]^2 - r[i]^2)
end

@kernel function _compute_spherical_centroids!(rc, rnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I

  r₀ = rnode[i]
  r₁ = rnode[i + 1]

  rc[i] = (3 / 4) * ((r₁^4 - r₀^4) / (r₁^3 - r₀^3))
end

@kernel function _compute_spherical_face_areas!(Aᵢ₊½, r, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, = I.I
  Aᵢ₊½[I] = 4π * r[i + 1]^2
end

@kernel function _compute_spherical_volumes!(V, r, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, = I.I
  V[I] = (4 / 3) * π * (r[i + 1]^3 - r[i]^3)
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
