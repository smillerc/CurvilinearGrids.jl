struct SphericalBasisCurvilinearGrid3D{CO,CNC,CE,NV,EM,CM,DL,CI,DS} <:
       AbstractCurvilinearGrid3D
  node_coordinates::CO
  cartesian_node_coordinates::CNC
  centroid_coordinates::CE
  node_velocities::NV
  edge_metrics::EM
  cell_center_metrics::CM
  nhalo::Int
  nnodes::NTuple{3,Int}
  domain_limits::DL
  iterators::CI
  discretization_scheme::DS
  is_static::Bool
  is_orthogonal::Bool
  discretization_scheme_name::Symbol
  halo_coords_included::Bool
end

function SphericalBasisCurvilinearGrid3D(
  r::AbstractArray{T,3},
  θ::AbstractArray{T,3},
  ϕ::AbstractArray{T,3},
  discretization_scheme::Symbol;
  backend=CPU(),
  is_static=false,
  is_orthogonal=false,
  init_metrics=true,
  empty_metrics=false,
  halo_coords_included=false,
  kwargs...,
) where {T}
  GradientDiscretizationScheme, order, use_symmetric_conservative_metric_scheme, nhalo, scheme_name = get_gradient_discretization_scheme(
    discretization_scheme
  )

  limits, iterators = get_iterators(size(r), halo_coords_included, nhalo)
  celldims = size(iterators.cell.full)

  discr_scheme = GradientDiscretizationScheme(
    order;
    use_cache=true,
    celldims=celldims,
    backend=backend,
    T=T,
    use_symmetric_conservative_metric_scheme=use_symmetric_conservative_metric_scheme,
  )

  coords, centroids, node_velocities, nnodes, cartesian_node_coords = _spherical_basis_grid_constructor(
    r, θ, ϕ, iterators, halo_coords_included; backend=backend
  )
  cell_center_metrics, edge_metrics = get_metric_soa(celldims, backend, T)

  m = SphericalBasisCurvilinearGrid3D(
    coords,
    cartesian_node_coords,
    centroids,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    nhalo,
    nnodes,
    limits,
    iterators,
    discr_scheme,
    is_static,
    is_orthogonal,
    scheme_name,
    halo_coords_included,
  )

  if init_metrics && !empty_metrics
    update!(m, true, halo_coords_included)
  end

  return m
end

function SphericalBasisCurvilinearGrid3D(
  r::AbstractVector{T},
  θ::AbstractVector{T},
  ϕ::AbstractVector{T},
  discretization_scheme::Symbol;
  kwargs...,
) where {T}
  ni = length(r)
  nj = length(θ)
  nk = length(ϕ)

  r3d = zeros(T, ni, nj, nk)
  θ3d = zeros(T, ni, nj, nk)
  ϕ3d = zeros(T, ni, nj, nk)

  @inbounds for k in 1:nk
    for j in 1:nj
      for i in 1:ni
        r3d[i, j, k] = r[i]
        θ3d[i, j, k] = θ[j]
        ϕ3d[i, j, k] = ϕ[k]
      end
    end
  end

  return SphericalBasisCurvilinearGrid3D(r3d, θ3d, ϕ3d, discretization_scheme; kwargs...)
end

function SphericalBasisCurvilinearGrid3D(
  (r0, θ0, ϕ0),
  (r1, θ1, ϕ1),
  (ni_cells, nj_cells, nk_cells)::NTuple{3,Int},
  discretization_scheme::Symbol;
  T=Float64,
  kwargs...,
)
  return SphericalBasisCurvilinearGrid3D(
    range(r0, r1; length=ni_cells + 1) .|> T,
    range(θ0, θ1; length=nj_cells + 1) .|> T,
    range(ϕ0, ϕ1; length=nk_cells + 1) .|> T,
    discretization_scheme;
    T=T,
    kwargs...,
  )
end

function _spherical_basis_grid_constructor(
  r::AbstractArray{T,3},
  θ::AbstractArray{T,3},
  ϕ::AbstractArray{T,3},
  domain_iterators,
  halo_coords_included;
  backend=KernelAbstractions.CPU(),
  kwargs...,
) where {T}
  celldims = size(domain_iterators.cell.full)
  nodedims = size(domain_iterators.node.full)

  coords = StructArray((
    r=KernelAbstractions.zeros(backend, T, nodedims),
    θ=KernelAbstractions.zeros(backend, T, nodedims),
    ϕ=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  centroids = StructArray((
    r=KernelAbstractions.zeros(backend, T, celldims),
    θ=KernelAbstractions.zeros(backend, T, celldims),
    ϕ=KernelAbstractions.zeros(backend, T, celldims),
  ))

  if halo_coords_included
    copy!(coords.r, r)
    copy!(coords.θ, θ)
    copy!(coords.ϕ, ϕ)
  else
    @views begin
      copy!(coords.r[domain_iterators.node.domain], r)
      copy!(coords.θ[domain_iterators.node.domain], θ)
      copy!(coords.ϕ[domain_iterators.node.domain], ϕ)
    end
  end

  if halo_coords_included
    domain = domain_iterators.cell.full
  else
    domain = domain_iterators.cell.domain
  end

  _spherical_basis_centroid_coordinates_kernel!(backend)(
    centroids, coords, domain; ndrange=size(domain)
  )

  node_velocities = StructArray((
    r=KernelAbstractions.zeros(backend, T, nodedims),
    θ=KernelAbstractions.zeros(backend, T, nodedims),
    ϕ=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  cartesian_node_coords = StructArray((
    x=KernelAbstractions.zeros(backend, T, nodedims),
    y=KernelAbstractions.zeros(backend, T, nodedims),
    z=KernelAbstractions.zeros(backend, T, nodedims),
  ))

  compute_xyz_coords!(
    cartesian_node_coords, coords, domain_iterators, backend, halo_coords_included
  )

  nnodes = size(r)
  return (coords, centroids, node_velocities, nnodes, cartesian_node_coords)
end

function update!(
  mesh::SphericalBasisCurvilinearGrid3D, force::Bool, include_halo_region::Bool
)
  if include_halo_region
    metric_domain = mesh.iterators.cell.full
  else
    metric_domain = mesh.iterators.cell.domain
  end

  if !mesh.is_static || force
    backend = KernelAbstractions.get_backend(mesh.centroid_coordinates.r)
    _spherical_basis_centroid_coordinates_kernel!(backend)(
      mesh.centroid_coordinates,
      mesh.node_coordinates,
      metric_domain;
      ndrange=size(metric_domain),
    )
    update_metrics!(mesh; include_halo_region=include_halo_region)
    _check_valid_metrics(mesh)
    compute_cartesian_node_coordinates!(mesh; include_halo_region=include_halo_region)
  else
    @warn("Attempting to update grid metrics when grid.is_static = true!")
  end
  return nothing
end

function compute_cartesian_node_coordinates!(
  mesh::SphericalBasisCurvilinearGrid3D; include_halo_region::Bool=mesh.halo_coords_included
)
  backend = KernelAbstractions.get_backend(mesh.node_coordinates.r)

  @kernel function _compute_xyz_coords!(x, y, z, r, θ, ϕ, domain)
    idx = @index(Global, Linear)
    I = domain[idx]

    rr = r[I]
    th = θ[I]
    ph = ϕ[I]

    x[I] = rr * sin(th) * cos(ph)
    y[I] = rr * sin(th) * sin(ph)
    z[I] = rr * cos(th)
  end

  dom = mesh.iterators.node.full

  _compute_xyz_coords!(backend)(
    mesh.cartesian_node_coordinates.x,
    mesh.cartesian_node_coordinates.y,
    mesh.cartesian_node_coordinates.z,
    mesh.node_coordinates.r,
    mesh.node_coordinates.θ,
    mesh.node_coordinates.ϕ,
    dom;
    ndrange=size(dom),
  )

  return nothing
end

function jacobian_matrix(mesh::SphericalBasisCurvilinearGrid3D, (i, j, k))
  xξ = mesh.cell_center_metrics.x₁.ξ
  yξ = mesh.cell_center_metrics.x₂.ξ
  zξ = mesh.cell_center_metrics.x₂.ξ
  xη = mesh.cell_center_metrics.x₁.η
  yη = mesh.cell_center_metrics.x₂.η
  zη = mesh.cell_center_metrics.x₂.η
  xζ = mesh.cell_center_metrics.x₁.η
  yζ = mesh.cell_center_metrics.x₂.η
  zζ = mesh.cell_center_metrics.x₂.η

  return @SMatrix [
    xξ[i, j, k] xη[i, j, k] xζ[i, j, k]
    yξ[i, j, k] yη[i, j, k] yζ[i, j, k]
    zξ[i, j, k] zη[i, j, k] zζ[i, j, k]
  ]
end

@kernel inbounds = true function _spherical_basis_centroid_coordinates_kernel!(
  centroids::StructArray{T,3}, coords::StructArray{T,3}, domain
) where {T}
  idx = @index(Global, Linear)
  didx = domain[idx]
  i, j, k = didx.I

  r = coords.r
  θ = coords.θ
  ϕ = coords.ϕ

  centroids.r[didx] =
    (
      r[i, j, k] +
      r[i + 1, j, k] +
      r[i + 1, j + 1, k] +
      r[i, j + 1, k] +
      r[i, j, k + 1] +
      r[i + 1, j, k + 1] +
      r[i + 1, j + 1, k + 1] +
      r[i, j + 1, k + 1]
    ) / 8

  centroids.θ[didx] =
    (
      θ[i, j, k] +
      θ[i + 1, j, k] +
      θ[i + 1, j + 1, k] +
      θ[i, j + 1, k] +
      θ[i, j, k + 1] +
      θ[i + 1, j, k + 1] +
      θ[i + 1, j + 1, k + 1] +
      θ[i, j + 1, k + 1]
    ) / 8

  centroids.ϕ[didx] =
    (
      ϕ[i, j, k] +
      ϕ[i + 1, j, k] +
      ϕ[i + 1, j + 1, k] +
      ϕ[i, j + 1, k] +
      ϕ[i, j, k + 1] +
      ϕ[i + 1, j, k + 1] +
      ϕ[i + 1, j + 1, k + 1] +
      ϕ[i, j + 1, k + 1]
    ) / 8
end
