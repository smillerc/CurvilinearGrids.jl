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
  r::AbstractVector{T}, θ::AbstractVector{T}, ϕ::AbstractVector{T}, nhalo::Int, backend
) where {T<:Real}

  #
  nr, ntheta, nphi = (length(r) + 2nhalo, length(θ) + 2nhalo, length(ϕ) + 2nhalo)
  halo_coords_included = true # the node dims include halo
  limits, domain_iterators = get_iterators((nr, ntheta, nphi), halo_coords_included, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  spherical_node_coords = (;
    r=KernelAbstractions.zeros(backend, T, nr),
    θ=KernelAbstractions.zeros(backend, T, ntheta),
    ϕ=KernelAbstractions.zeros(backend, T, nphi),
  )

  @views begin
    copy!(spherical_node_coords.r[domain_iterators.node.domain.indices[1]], r)
    copy!(spherical_node_coords.θ[domain_iterators.node.domain.indices[2]], θ)
    copy!(spherical_node_coords.ϕ[domain_iterators.node.domain.indices[3]], ϕ)
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
    cartesian_node_coords, spherical_node_coords, domain_iterators, backend
  )
  compute_centroids!(
    spherical_centroid_coords, spherical_node_coords, domain_iterators, backend
  )
  compute_volumes!(cell_volumes, spherical_node_coords, domain_iterators, backend)

  compute_face_areas!(face_areas, spherical_node_coords, domain_iterators, backend)

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

function compute_centroids!(centroids, node_coordinates, iters, backend)
  kern = _compute_centroids!(backend)
  kern(
    centroids.r,
    centroids.θ,
    centroids.ϕ,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    iters.cell.domain;
    ndrange=size(iters.cell.domain),
  )
  return nothing
end

function compute_face_areas!(face_areas, node_coordinates, iters, backend)
  kern = _compute_face_areas!(backend)
  kern(
    face_areas.i₊½,
    face_areas.j₊½,
    face_areas.k₊½,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    iters.cell.domain;
    ndrange=size(iters.cell.domain),
  )
  return nothing
end

function compute_volumes!(volumes, node_coordinates, iters, backend)
  kern = _compute_volumes!(backend)
  kern(
    volumes,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    iters.cell.domain;
    ndrange=size(iters.cell.domain),
  )
  return nothing
end

function compute_xyz_coords!(xyz_coordinates, node_coordinates, iters, backend)
  kern = _compute_xyz_coords!(backend)
  kern(
    xyz_coordinates.x,
    xyz_coordinates.y,
    xyz_coordinates.z,
    node_coordinates.r,
    node_coordinates.θ,
    node_coordinates.ϕ,
    iters.node.domain;
    ndrange=size(iters.node.domain),
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

@kernel function _compute_face_areas!(Aᵢ₊½, Aⱼ₊½, Aₖ₊½, r, θ, ϕ, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j, k = I.I

  Δμ = cos(θ[j]) - cos(θ[j + 1])
  Δϕ = ϕ[k + 1] - ϕ[k]
  Δr² = r[i + 1]^2 - r[i]^2

  Aᵢ₊½[I] = r[i]^2 * Δμ * Δϕ

  Aⱼ₊½[I] = (1 / 2) * Δr² * sin(θ[j]) * Δϕ
  Aₖ₊½[I] = (1 / 2) * Δr² * Δμ
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
