function SphericalGrid3D(
  _r::AbstractVector{T},
  _θ::AbstractVector{T},
  _ϕ::AbstractVector{T},
  nhalo::Int,
  backend=CPU();
  halo_coords_included=false,
) where {T<:Real}

  # This pads the coordinate arrays with halo geometry using the 
  # spacing on either end of the coordinate vector. For this mesh type, 
  # halo geometry MUST be defined
  if !halo_coords_included
    r = pad_with_halo(_r, nhalo)
    θ = pad_with_halo(_θ, nhalo)
    ϕ = pad_with_halo(_ϕ, nhalo)
    halo_coords_included = true
  else
    r = _r
    θ = _θ
    ϕ = _ϕ
  end

  nr, ntheta, nphi = length(r), length(θ), length(ϕ)

  limits, domain_iterators = get_iterators((nr, ntheta, nphi), halo_coords_included, nhalo)

  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  spherical_node_coords = (
    KernelAbstractions.zeros(backend, T, nodedims[1]),
    KernelAbstractions.zeros(backend, T, nodedims[2]),
    KernelAbstractions.zeros(backend, T, nodedims[3]),
  )

  if halo_coords_included
    @views begin
      copy!(spherical_node_coords[1], r)
      copy!(spherical_node_coords[2], θ)
      copy!(spherical_node_coords[3], ϕ)
    end
  else
    @views begin
      copy!(spherical_node_coords[1][domain_iterators.node.domain.indices[1]], r)
      copy!(spherical_node_coords[2][domain_iterators.node.domain.indices[2]], θ)
      copy!(spherical_node_coords[3][domain_iterators.node.domain.indices[3]], ϕ)
    end
  end

  spherical_centroid_coords = (
    KernelAbstractions.zeros(backend, T, celldims[1]),
    KernelAbstractions.zeros(backend, T, celldims[2]),
    KernelAbstractions.zeros(backend, T, celldims[3]),
  )

  cell_volumes = KernelAbstractions.zeros(backend, T, celldims)

  face_areas = (
    KernelAbstractions.zeros(backend, T, celldims),
    KernelAbstractions.zeros(backend, T, celldims),
    KernelAbstractions.zeros(backend, T, celldims),
  )

  mesh = OrthogonalGrid{
    3,
    T,
    SphericalCS,
    typeof(spherical_node_coords),
    typeof(spherical_centroid_coords),
    typeof(cell_volumes),
    typeof(domain_iterators),
    typeof(limits),
    typeof(face_areas),
  }(
    spherical_node_coords,
    spherical_centroid_coords,
    cell_volumes,
    domain_iterators,
    limits,
    face_areas,
    nhalo,
  )

  update!(mesh)

  return mesh
end

function update!(mesh::OrthogonalGrid{3,T,SphericalCS}) where {T}
  backend = KernelAbstractions.get_backend(mesh.node_coordinates[1])
  halo_coords_included = true

  compute_centroids!(
    mesh.centroid_coordinates,
    mesh.node_coordinates,
    mesh.iterators,
    backend,
    halo_coords_included,
  )

  compute_volumes!(
    mesh.cell_volumes, mesh.node_coordinates, mesh.iterators, backend, halo_coords_included
  )

  compute_face_areas!(
    mesh.face_areas,
    mesh.node_coordinates,
    mesh.iterators,
    backend,
    halo_coords_included,
    mesh.nhalo,
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

function compute_face_areas!(
  face_areas, node_coordinates, iters, backend, halo_coords_included, nhalo
)
  domain = iters.cell.domain
  if halo_coords_included && nhalo > 0
    offset = +nhalo
    i₊½_domain = expand_upper(expand_lower(domain, 1, offset), 1, offset - 1)
    j₊½_domain = expand_upper(expand_lower(domain, 2, offset), 2, offset - 1)
    k₊½_domain = expand_upper(expand_lower(domain, 3, offset), 3, offset - 1)
  elseif halo_coords_included
    i₊½_domain = domain
    j₊½_domain = domain
    k₊½_domain = domain
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
    face_areas[1],
    node_coordinates[1],
    node_coordinates[2],
    node_coordinates[3],
    i₊½_domain;
    ndrange=size(i₊½_domain),
  )

  theta_kernel(
    face_areas[2],
    node_coordinates[1],
    node_coordinates[2],
    node_coordinates[3],
    j₊½_domain;
    ndrange=size(j₊½_domain),
  )

  phi_kernel(
    face_areas[3],
    node_coordinates[1],
    node_coordinates[2],
    node_coordinates[3],
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
    node_coordinates[1],
    node_coordinates[2],
    node_coordinates[3],
    domain;
    ndrange=size(domain),
  )
  return nothing
end

@kernel function _compute_centroids!(rc, θc, ϕc, rnode, θnode, ϕnode, cell_domain)
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j, k = I.I

  r1 = rnode[i]
  r2 = rnode[i + 1]
  θ₀ = θnode[j]
  θ₁ = θnode[j + 1]
  ϕ₀ = ϕnode[k]
  ϕ₁ = ϕnode[k + 1]

  num = r2^3 + r2^2 * r1 + r2 * r1^2 + r1^3
  den = r2^2 + r2 * r1 + r1^2
  rc[i] = (3 / 4) * (num / den)

  # rc[i] = (3 / 4) * ((r2^4 - r1^4) / (r2^3 - r1^3))
  θc[j] = acos((cos(θ₀) + cos(θ₁)) / 2)
  ϕc[k] = (ϕ₀ + ϕ₁) / 2
end

@kernel function _compute_radial_face_areas!(
  Aᵢ₊½::AbstractArray{T,N}, r, θ, ϕ, domain
) where {T,N}
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I # these are CELL indices

  # these are NODE indexed
  Δμ = cos(θ[j]) - cos(θ[j + 1])
  Δϕ = ϕ[k + 1] - ϕ[k]

  # this uses CELL indexing, e.g. for cell i,j,k, the i+1/2 face area is...
  Aᵢ₊½[i, j, k] = r[i + 1]^2 * Δμ * Δϕ
end

@kernel function _compute_theta_face_areas!(
  Aⱼ₊½::AbstractArray{T,N}, r, θ, ϕ, domain
) where {T,N}
  idx = @index(Global, Linear)
  I = domain[idx]
  i, j, k = I.I # these are CELL indices

  # these are NODE indexed
  Δϕ = ϕ[k + 1] - ϕ[k]
  Δr² = r[i + 1]^2 - r[i]^2

  # this uses CELL indexing, e.g. for cell i,j,k, the j+1/2 face area is...
  Aⱼ₊½[i, j, k] = T(1 / 2) * Δr² * sin(θ[j + 1]) * Δϕ
end

@kernel function _compute_phi_face_areas!(
  Aₖ₊½::AbstractArray{T,N}, r, θ, ϕ, domain
) where {T,N}
  idx = @index(Global, Linear)
  I = domain[idx] # these are CELL indices
  i, j, k = I.I

  # these are NODE indexed
  Δμ = cos(θ[j]) - cos(θ[j + 1])
  Δr² = r[i + 1]^2 - r[i]^2

  # this uses CELL indexing, e.g. for cell i,j,k, the k+1/2 face area is...
  Aₖ₊½[i, j, k] = T(1 / 2) * Δr² * Δμ
end

@kernel function _compute_volumes!(V::AbstractArray{T,N}, r, θ, ϕ, cell_domain) where {T,N}
  idx = @index(Global, Linear)
  I = cell_domain[idx]
  i, j, k = I.I

  # r2 = r[i + 1]
  # r1 = r[i]
  # Δr³ = ((r2 - r1) * (r2^2 + r2 * r1 + r1^2)) / 3
  # # Δr³ = (r[i + 1]^3 - r[i]^3) / 3
  # Δμ = cos(θ[j]) - cos(θ[j + 1])
  # Δϕ = ϕ[k + 1] - ϕ[k]
  # V[I] = Δr³ * Δμ * Δϕ

  Fr = radial_factor_stable(r[i], r[i + 1])
  Fθ = theta_factor_stable(θ[j], θ[j + 1])
  Fφ = ϕ[k + 1] - ϕ[k]

  V[I] = Fr * Fθ * Fφ
end

@inline function sinx_over_x(x::T) where {T<:AbstractFloat}
  ax = abs(x)
  if ax < T(1e-6)
    # 1 - x^2/6 + x^4/120 is plenty for Float64/Float32 stability here
    x2 = x * x
    return one(T) - x2 / T(6) + (x2 * x2) / T(120)
  else
    return sin(x) / x
  end
end

@inline function dphi_0_2pi(φ1::T, φ2::T) where {T<:AbstractFloat}
  Δ = φ2 - φ1
  return (Δ >= zero(T)) ? Δ : (Δ + T(2pi))
end

@inline function radial_factor_stable(r1::T, r2::T) where {T<:AbstractFloat}
  Δr = r2 - r1
  # s = r2^2 + r2*r1 + r1^2 with a couple fma/muladd steps for better rounding
  s = muladd(r2, r2, r1 * r1)   # r2^2 + r1^2
  s = muladd(r2, r1, s)       # + r2*r1
  return (Δr * s) / T(3)       # (r2^3 - r1^3)/3
end

@inline function theta_factor_stable(θ1::T, θ2::T) where {T<:AbstractFloat}
  θm = (θ1 + θ2) / T(2)
  h = (θ2 - θ1) / T(2)
  # cosθ1 - cosθ2 = 2*sin(θm)*sin(h) = 2*sin(θm)*(h*sinc(h))
  return (T(2) * sin(θm)) * (h * sinx_over_x(h))
end

function face_location(
  mesh::OrthogonalGrid{3,T,SphericalCS}, I::CartesianIndex{3}, axis
) where {T}
  i, j, k = I.I
  if axis == 1
    return @SVector [
      mesh.centroid_coordinates[1][i + 1],
      mesh.node_coordinates[2][j],
      mesh.centroid_coordinates[3][k],
    ]
  elseif axis == 2
    return @SVector [
      mesh.centroid_coordinates[1][i],
      mesh.node_coordinates[2][j + 1],
      mesh.centroid_coordinates[3][k],
    ]
  else
    return @SVector [
      mesh.centroid_coordinates[1][i],
      mesh.centroid_coordinates[2][j],
      mesh.node_coordinates[3][k + 1],
    ]
  end
end
