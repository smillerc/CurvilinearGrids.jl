
struct SphericalGrid3D{NC,CNC,CC,CV,I,DL,FA} <: AbstractCurvilinearGrid3D
  node_coordinates::NC
  cartesian_node_coordinates::CNC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

function SphericalGrid3D(
  r::AbstractVector{T},
  θ::AbstractVector{T},
  ϕ::AbstractVector{T},
  nhalo::Int,
  backend;
  halo_coords_included=false,
) where {T<:Real}
  nr, ntheta, nphi = length(r), length(θ), length(ϕ)
  limits, domain_iterators = get_iterators((nr, ntheta, nphi), halo_coords_included, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  spherical_node_coords = (;
    r=KernelAbstractions.zeros(backend, T, nodedims[1]),
    θ=KernelAbstractions.zeros(backend, T, nodedims[2]),
    ϕ=KernelAbstractions.zeros(backend, T, nodedims[3]),
  )

  if halo_coords_included
    @views begin
      copy!(spherical_node_coords.r, r)
      copy!(spherical_node_coords.θ, θ)
      copy!(spherical_node_coords.ϕ, ϕ)
    end
  else
    @views begin
      copy!(spherical_node_coords.r[domain_iterators.node.domain.indices[1]], r)
      copy!(spherical_node_coords.θ[domain_iterators.node.domain.indices[2]], θ)
      copy!(spherical_node_coords.ϕ[domain_iterators.node.domain.indices[3]], ϕ)
    end
  end

  cartesian_node_coords = (;
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  )

  spherical_centroid_coords = (;
    r=KernelAbstractions.zeros(backend, T, celldims[1]),
    θ=KernelAbstractions.zeros(backend, T, celldims[2]),
    ϕ=KernelAbstractions.zeros(backend, T, celldims[3]),
  )

  cell_volumes = KernelAbstractions.zeros(backend, T, celldims)

  face_areas = (;
    i₊½=KernelAbstractions.zeros(backend, T, celldims),
    j₊½=KernelAbstractions.zeros(backend, T, celldims),
    k₊½=KernelAbstractions.zeros(backend, T, celldims),
  )

  compute_xyz_coords!(
    cartesian_node_coords,
    spherical_node_coords,
    domain_iterators,
    backend,
    halo_coords_included,
  )
  compute_centroids!(
    spherical_centroid_coords,
    spherical_node_coords,
    domain_iterators,
    backend,
    halo_coords_included,
  )
  compute_volumes!(
    cell_volumes, spherical_node_coords, domain_iterators, backend, halo_coords_included
  )

  compute_face_areas!(
    face_areas,
    spherical_node_coords,
    domain_iterators,
    backend,
    halo_coords_included,
    nhalo,
  )

  return SphericalGrid3D(
    spherical_node_coords,
    cartesian_node_coords,
    spherical_centroid_coords,
    cell_volumes,
    domain_iterators,
    limits,
    face_areas,
    nhalo,
  )
end

function compute_centroids!(
  centroids, node_coordinates, iters, backend, halo_coords_included
)
  if halo_coords_included
    domain = iters.cell.full
  else
    domain = iters.cell.domain
  end

  kern = _compute_centroids!(backend)
  kern(
    centroids.r,
    centroids.θ,
    centroids.ϕ,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    domain;
    ndrange=size(domain),
  )
  return nothing
end

function compute_face_areas!(
  face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
)
  domain = iters.cell.domain
  if halo_coords_included
    offset = +nhalo
    i₊½_domain = expand_upper(expand_lower(domain, 1, offset), 1, offset - 1)
    j₊½_domain = expand_upper(expand_lower(domain, 2, offset), 2, offset - 1)
    k₊½_domain = expand_upper(expand_lower(domain, 3, offset), 3, offset - 1)
  else
    offset = +1
    i₊½_domain = expand_lower(domain, 1, offset)
    j₊½_domain = expand_lower(domain, 2, offset)
    k₊½_domain = expand_lower(domain, 3, offset)
  end

  r_kernel = _compute_radial_face_areas!(backend)
  theta_kernel = _compute_theta_face_areas!(backend)
  phi_kernel = _compute_phi_face_areas!(backend)

  r_kernel(
    face_areas.i₊½,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    i₊½_domain;
    ndrange=size(i₊½_domain),
  )

  theta_kernel(
    face_areas.j₊½,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    j₊½_domain;
    ndrange=size(j₊½_domain),
  )

  phi_kernel(
    face_areas.k₊½,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    k₊½_domain;
    ndrange=size(k₊½_domain),
  )
  return nothing
end

function compute_volumes!(volumes, node_coordinates, iters, backend, halo_coords_included)
  if halo_coords_included
    domain = iters.cell.full
  else
    domain = iters.cell.domain
  end

  kern = _compute_volumes!(backend)
  kern(
    volumes,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    domain;
    ndrange=size(domain),
  )
  return nothing
end

function compute_xyz_coords!(
  xyz_coordinates, node_coordinates, iters, backend, halo_coords_included
)
  if halo_coords_included
    domain = iters.node.full
  else
    domain = iters.node.domain
  end

  kern = _compute_xyz_coords!(backend)
  kern(
    xyz_coordinates.x,
    xyz_coordinates.y,
    xyz_coordinates.z,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    domain;
    ndrange=size(domain),
  )
  return nothing
end

@kernel function _compute_centroids!(rc, θc, ϕc, rnode, θnode, ϕnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j, k = I.I

  r₀ = rnode[i]
  r₁ = rnode[i + 1]
  θ₀ = θnode[j]
  θ₁ = θnode[j + 1]
  ϕ₀ = ϕnode[k]
  ϕ₁ = ϕnode[k + 1]

  rc[i] = (3 / 4) * ((r₁^4 - r₀^4) / (r₁^3 - r₀^3))
  θc[j] = acos((cos(θ₀) + cos(θ₁)) / 2)
  ϕc[k] = (ϕ₀ + ϕ₁) / 2
end

@kernel function _compute_radial_face_areas!(Aᵢ₊½, r, θ, ϕ, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I # these are CELL indices

  # these are NODE indexed
  Δμ = cos(θ[j]) - cos(θ[j + 1])
  Δϕ = ϕ[k + 1] - ϕ[k]

  # this uses CELL indexing, e.g. for cell i,j,k, the i+1/2 face area is...
  Aᵢ₊½[i, j, k] = r[i + 1]^2 * Δμ * Δϕ
end

@kernel function _compute_theta_face_areas!(Aⱼ₊½, r, θ, ϕ, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I # these are CELL indices

  # these are NODE indexed
  Δϕ = ϕ[k + 1] - ϕ[k]
  Δr² = r[i + 1]^2 - r[i]^2

  # this uses CELL indexing, e.g. for cell i,j,k, the j+1/2 face area is...
  Aⱼ₊½[i, j, k] = (1 / 2) * Δr² * sin(θ[j + 1]) * Δϕ
end

@kernel function _compute_phi_face_areas!(Aₖ₊½, r, θ, ϕ, domain)
  idx = @index(Global, Linear)
  I = domain[idx] # these are CELL indices
  i, j, k = I.I

  # these are NODE indexed
  Δμ = cos(θ[j]) - cos(θ[j + 1])
  Δr² = r[i + 1]^2 - r[i]^2

  # this uses CELL indexing, e.g. for cell i,j,k, the k+1/2 face area is...
  Aₖ₊½[i, j, k] = (1 / 2) * Δr² * Δμ
end

@kernel function _compute_volumes!(V, r, θ, ϕ, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j, k = I.I

  Δr³ = r[i + 1]^3 - r[i]^3
  Δμ = cos(θ[j]) - cos(θ[j + 1])
  Δϕ = ϕ[k + 1] - ϕ[k]

  V[I] = (1 / 3) * Δr³ * Δμ * Δϕ
end

@kernel function _compute_xyz_coords!(x, y, z, r, θ, ϕ, domain)
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I

  rr = r[i]
  th = θ[j]
  ph = ϕ[k]

  x[I] = rr * sin(th) * cos(ph)
  y[I] = rr * sin(th) * sin(ph)
  z[I] = rr * cos(th)
end

function face_location(mesh::SphericalGrid3D, I::CartesianIndex{3}, axis)
  i, j, k = I.I
  if axis == 1
    return @SVector [
      mesh.centroid_coordinates.r[i + 1],
      mesh.node_coordinates.θ[j],
      mesh.centroid_coordinates.ϕ[k],
    ]
  elseif axis == 2
    return @SVector [
      mesh.centroid_coordinates.r[i],
      mesh.node_coordinates.θ[j + 1],
      mesh.centroid_coordinates.ϕ[k],
    ]
  else
    return @SVector [
      mesh.centroid_coordinates.r[i],
      mesh.centroid_coordinates.θ[j],
      mesh.node_coordinates.ϕ[k + 1],
    ]
  end
end
