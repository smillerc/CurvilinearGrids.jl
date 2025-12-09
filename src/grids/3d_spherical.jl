
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

# --------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------

@inline cell_volume(mesh::SphericalGrid3D, I::CartesianIndex{3}) =
  @inbounds mesh.cell_volumes[I]

@inline function face_area_p(mesh::SphericalGrid3D, I::CartesianIndex{3}, axis::Int)
  @inbounds if axis == 1
    return mesh.face_areas.i₊½[I]
  elseif axis == 2
    return mesh.face_areas.j₊½[I]
  else
    return mesh.face_areas.k₊½[I]
  end
end

@inline function face_area_m(mesh::SphericalGrid3D, I::CartesianIndex{3}, axis::Int)
  I₋ = CartesianDomains.shift(I, axis, -1)
  @inbounds if axis == 1
    return mesh.face_areas.i₊½[I₋]  # i-1/2
  elseif axis == 2
    return mesh.face_areas.j₊½[I₋]  # j-1/2
  else
    return mesh.face_areas.k₊½[I₋]  # k-1/2
  end
end

@inline function face_val(A, I₁::CartesianIndex{3}, I₂::CartesianIndex{3})
  @inbounds (1 / 2) * (A[I₁] + A[I₂])
end

# Effective physical spacing at edge between I and I+ê_axis
@inline function Δx_eff_edge(mesh::SphericalGrid3D, I::CartesianIndex{3}, axis::Int)
  I₊ = CartesianDomains.shift(I, axis, +1)
  V_L = cell_volume(mesh, I)
  V_R = cell_volume(mesh, I₊)
  A_n = face_area_p(mesh, I, axis)
  return (V_L + V_R) / A_n
end

# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

function cell_center_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}
  I₊ = CartesianDomains.shift(I, axis, +1)
  I₋ = CartesianDomains.shift(I, axis, -1)

  @inbounds A₊ = A[I₊]
  @inbounds Aᵢ = A[I]
  @inbounds A₋ = A[I₋]

  Vᵢ = cell_volume(mesh, I)
  V₊ = cell_volume(mesh, I₊)
  V₋ = cell_volume(mesh, I₋)

  A_face_p = face_area_p(mesh, I, axis)
  A_face_m = face_area_m(mesh, I, axis)

  Δx_f = (Vᵢ + V₊) / A_face_p    # forward spacing
  Δx_b = (V₋ + Vᵢ) / A_face_m    # backward spacing

  ∂A_f = (A₊ - Aᵢ) / Δx_f
  ∂A_b = (Aᵢ - A₋) / Δx_b

  ∂A = (1 / 2) * (∂A_f + ∂A_b)
  return ∂A
end

# ∇A at the cell center (physical gradient)
function cell_center_gradient(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∂A∂r = cell_center_derivative(mesh, A, I, 1)
  ∂A∂θ = cell_center_derivative(mesh, A, I, 2)
  ∂A∂ϕ = cell_center_derivative(mesh, A, I, 3)
  return @SVector [∂A∂r, ∂A∂θ, ∂A∂ϕ]
end

# ∇⋅A at the cell center (A_r,A_θ,A_ϕ are physical components)
function cell_center_divergence(
  mesh::SphericalGrid3D,
  (A_r::AbstractArray{T,3}, A_θ::AbstractArray{T,3}, A_ϕ::AbstractArray{T,3}),
  I::CartesianIndex{3},
) where {T}
  V = cell_volume(mesh, I)

  # convenience neighbor indices
  I_i₊ = CartesianDomains.shift(I, 1, +1)
  I_i₋ = CartesianDomains.shift(I, 1, -1)
  I_j₊ = CartesianDomains.shift(I, 2, +1)
  I_j₋ = CartesianDomains.shift(I, 2, -1)
  I_k₊ = CartesianDomains.shift(I, 3, +1)
  I_k₋ = CartesianDomains.shift(I, 3, -1)

  # radial faces
  @inbounds A_rᵢ = A_r[I]
  @inbounds A_rᵢ₊₁ = A_r[I_i₊]
  @inbounds A_rᵢ₋₁ = A_r[I_i₋]

  A_r_face_p = (1 / 2) * (A_rᵢ + A_rᵢ₊₁)
  A_r_face_m = (1 / 2) * (A_rᵢ₋₁ + A_rᵢ)

  A_r_area_p = face_area_p(mesh, I, 1)
  A_r_area_m = face_area_m(mesh, I, 1)

  flux_r = A_r_face_p * A_r_area_p - A_r_face_m * A_r_area_m

  # theta faces
  @inbounds A_θⱼ = A_θ[I]
  @inbounds A_θⱼ₊₁ = A_θ[I_j₊]
  @inbounds A_θⱼ₋₁ = A_θ[I_j₋]

  A_θ_face_p = (1 / 2) * (A_θⱼ + A_θⱼ₊₁)
  A_θ_face_m = (1 / 2) * (A_θⱼ₋₁ + A_θⱼ)

  A_θ_area_p = face_area_p(mesh, I, 2)
  A_θ_area_m = face_area_m(mesh, I, 2)

  flux_θ = A_θ_face_p * A_θ_area_p - A_θ_face_m * A_θ_area_m

  # phi faces
  @inbounds A_ϕₖ = A_ϕ[I]
  @inbounds A_ϕₖ₊₁ = A_ϕ[I_k₊]
  @inbounds A_ϕₖ₋₁ = A_ϕ[I_k₋]

  A_ϕ_face_p = (1 / 2) * (A_ϕₖ + A_ϕₖ₊₁)
  A_ϕ_face_m = (1 / 2) * (A_ϕₖ₋₁ + A_ϕₖ)

  A_ϕ_area_p = face_area_p(mesh, I, 3)
  A_ϕ_area_m = face_area_m(mesh, I, 3)

  flux_ϕ = A_ϕ_face_p * A_ϕ_area_p - A_ϕ_face_m * A_ϕ_area_m

  div_A = (flux_r + flux_θ + flux_ϕ) / V
  return div_A
end

# ∇×A at the cell center (physical components)
function cell_center_curl(
  mesh::SphericalGrid3D,
  (A_r::AbstractArray{T,3}, A_θ::AbstractArray{T,3}, A_ϕ::AbstractArray{T,3}),
  I::CartesianIndex{3},
) where {T}
  V = cell_volume(mesh, I)

  I_i₊ = CartesianDomains.shift(I, 1, +1)
  I_i₋ = CartesianDomains.shift(I, 1, -1)
  I_j₊ = CartesianDomains.shift(I, 2, +1)
  I_j₋ = CartesianDomains.shift(I, 2, -1)
  I_k₊ = CartesianDomains.shift(I, 3, +1)
  I_k₋ = CartesianDomains.shift(I, 3, -1)

  # ω_r: circulation in (θ,ϕ)
  A_ϕ_face_θp = face_val(A_ϕ, I, I_j₊)
  A_ϕ_face_θm = face_val(A_ϕ, I_j₋, I)

  A_θ_face_ϕp = face_val(A_θ, I, I_k₊)
  A_θ_face_ϕm = face_val(A_θ, I_k₋, I)

  A_θ_area_p = face_area_p(mesh, I, 2)
  A_θ_area_m = face_area_m(mesh, I, 2)
  A_ϕ_area_p = face_area_p(mesh, I, 3)
  A_ϕ_area_m = face_area_m(mesh, I, 3)

  circ_r =
    (A_θ_area_p * A_ϕ_face_θp - A_θ_area_m * A_ϕ_face_θm) -
    (A_ϕ_area_p * A_θ_face_ϕp - A_ϕ_area_m * A_θ_face_ϕm)

  ω_r = circ_r / V

  # ω_θ: circulation in (ϕ,r)
  A_r_face_ϕp = face_val(A_r, I, I_k₊)
  A_r_face_ϕm = face_val(A_r, I_k₋, I)

  A_ϕ_face_rp = face_val(A_ϕ, I, I_i₊)
  A_ϕ_face_rm = face_val(A_ϕ, I_i₋, I)

  A_ϕ_area_p2 = face_area_p(mesh, I, 3)
  A_ϕ_area_m2 = face_area_m(mesh, I, 3)
  A_r_area_p = face_area_p(mesh, I, 1)
  A_r_area_m = face_area_m(mesh, I, 1)

  circ_θ =
    (A_ϕ_area_p2 * A_r_face_ϕp - A_ϕ_area_m2 * A_r_face_ϕm) -
    (A_r_area_p * A_ϕ_face_rp - A_r_area_m * A_ϕ_face_rm)

  ω_θ = circ_θ / V

  # ω_ϕ: circulation in (r,θ)
  A_θ_face_rp2 = face_val(A_θ, I, I_i₊)
  A_θ_face_rm2 = face_val(A_θ, I_i₋, I)

  A_r_face_θp2 = face_val(A_r, I, I_j₊)
  A_r_face_θm2 = face_val(A_r, I_j₋, I)

  A_r_area_p2 = face_area_p(mesh, I, 1)
  A_r_area_m2 = face_area_m(mesh, I, 1)
  A_θ_area_p2 = face_area_p(mesh, I, 2)
  A_θ_area_m2 = face_area_m(mesh, I, 2)

  circ_ϕ =
    (A_r_area_p2 * A_θ_face_rp2 - A_r_area_m2 * A_θ_face_rm2) -
    (A_θ_area_p2 * A_r_face_θp2 - A_θ_area_m2 * A_r_face_θm2)

  ω_ϕ = circ_ϕ / V

  return @SVector [ω_r, ω_θ, ω_ϕ]
end

# --------------------------------------------------------------------
# Edge operators
# --------------------------------------------------------------------

function edge_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}
  I₊ = CartesianDomains.shift(I, axis, +1)
  @inbounds Aᵢ = A[I]
  @inbounds A₊ = A[I₊]

  Δx = Δx_eff_edge(mesh, I, axis)
  ∂A = (A₊ - Aᵢ) / Δx
  return ∂A
end

# ∇A at the edge center (full physical gradient)
function edge_gradient(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∂A∂r = edge_derivative(mesh, A, I, 1)
  ∂A∂θ = edge_derivative(mesh, A, I, 2)
  ∂A∂ϕ = edge_derivative(mesh, A, I, 3)
  return @SVector [∂A∂r, ∂A∂θ, ∂A∂ϕ]
end

# ∇⋅A at the edge center = average of adjacent cell-center divergences
function edge_divergence(
  mesh::SphericalGrid3D,
  (A_r::AbstractArray{T,3}, A_θ::AbstractArray{T,3}, A_ϕ::AbstractArray{T,3}),
  I::CartesianIndex{3},
  edge_axis::Int,
) where {T}
  I₊ = CartesianDomains.shift(I, edge_axis, +1)

  div_L = cell_center_divergence(mesh, (A_r, A_θ, A_ϕ), I)
  div_R = cell_center_divergence(mesh, (A_r, A_θ, A_ϕ), I₊)

  return (1 / 2) * (div_L + div_R)
end

# ∇×A at the edge center = average of adjacent cell-center curls
function edge_curl(
  mesh::SphericalGrid3D,
  (
    A_r_edge::AbstractArray{T,3}, A_θ_edge::AbstractArray{T,3}, A_ϕ_edge::AbstractArray{T,3}
  ),
  I::CartesianIndex{3},
  edge_axis::Int,
) where {T}
  I₊ = CartesianDomains.shift(I, edge_axis, +1)

  curl_L = cell_center_curl(mesh, (A_r_edge, A_θ_edge, A_ϕ_edge), I)
  curl_R = cell_center_curl(mesh, (A_r_edge, A_θ_edge, A_ϕ_edge), I₊)

  return (1 / 2) * (curl_L + curl_R)
end
