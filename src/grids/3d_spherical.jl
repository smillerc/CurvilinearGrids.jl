module SphericalGrid

using CartesianDomains

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
  r_kernel = _compute_radial_face_areas!(backend)
  theta_kernel = _compute_theta_face_areas!(backend)
  phi_kernel = _compute_phi_face_areas!(backend)

  i₊½_domain = expand_lower(iters.cell.domain, 1, +1)
  j₊½_domain = expand_lower(iters.cell.domain, 2, +1)
  k₊½_domain = expand_lower(iters.cell.domain, 3, +1)

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

#
# Cell-center operators
#

function cell_center_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}

  #
  ᵢ₊₁ = CartesianDomains.shift(I, axis, +1)
  @inbounds Aᵢ = A[I]
  @inbounds Aᵢ₊₁ = A[ᵢ₊₁]
  # ∂A = ...
  return ∂A
end

# ∇A at the cell center
function cell_center_gradient(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∂A∂r = cell_center_derivative(mesh, A, I, 1)
  ∂A∂θ = cell_center_derivative(mesh, A, I, 2)
  ∂A∂ϕ = cell_center_derivative(mesh, A, I, 3)
  return @SVector [∂A∂r, ∂A∂θ, ∂A∂ϕ]
end

# ∇×A at the cell center
function cell_center_curl(
  mesh::SphericalGrid3D,
  (Ar::AbstractArray{T,3}, Aθ::AbstractArray{T,3}, Aϕ::AbstractArray{T,3}),
  I::CartesianIndex{3},
  axis::Int,
) where {T}
  ∇Ar = cell_center_gradient(mesh, Ar, I)
  ∇Aθ = cell_center_gradient(mesh, Aθ, I)
  ∇Aϕ = cell_center_gradient(mesh, Aϕ, I)
  return curl_A
end

# ∇⋅A at the cell center
function cell_center_divergence(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∇A = cell_center_gradient(mesh, A, I)
  return div_A
end

#
# Edge operators
#

function edge_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}

  #
  ᵢ₊₁ = CartesianDomains.shift(I, axis, +1)
  @inbounds Aᵢ = A[I]
  @inbounds Aᵢ₊₁ = A[ᵢ₊₁]
  # ∂A = ...
  return ∂A
end

# ∇A at the edge center
function edge_gradient(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∂A∂r = edge_derivative(mesh, A, I, 1)
  ∂A∂θ = edge_derivative(mesh, A, I, 2)
  ∂A∂ϕ = edge_derivative(mesh, A, I, 3)
  return @SVector [∂A∂r, ∂A∂θ, ∂A∂ϕ]
end

# ∇×A at the edge center
function edge_curl(
  mesh::SphericalGrid3D,
  (Ar_edge::AbstractArray{T,3}, Aθ_edge::AbstractArray{T,3}, Aϕ_edge::AbstractArray{T,3}),
  I::CartesianIndex{3},
  edge_axis::Int,
) where {T}
  ∇Ar_edge = edge_gradient(mesh, Ar_edge, I)
  ∇Aθ_edge = edge_gradient(mesh, Aθ_edge, I)
  ∇Aϕ_edge = edge_gradient(mesh, Aϕ_edge, I)
  return curl_A
end

# ∇⋅A at the cell center
function edge_divergence(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, edge_axis::Int
) where {T}
  ∇A = edge_gradient(mesh, A, I)
  return div_A
end

end # module