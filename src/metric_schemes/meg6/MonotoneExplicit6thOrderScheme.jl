module MonotoneExplicit6thOrderScheme

using LinearAlgebra
using StaticArrays, MappedArrays
using KernelAbstractions

using ..IndexingUtils

export MonotoneExplicit6thOrderDiscretization
export update_edge_conserved_metrics!, update_cell_center_metrics!

const ξ = 1
const η = 2
const ζ = 3

include("gradients.jl")
include("interpolation.jl")

struct MonotoneExplicit6thOrderDiscretization{C}
  cache::C
end

function MonotoneExplicit6thOrderDiscretization(
  domain::CartesianIndices{N}, T=Float64
) where {N}
  cellsize = size(domain)

  if any(cellsize .<= 3)
    error(
      "Domain size ($(cellsize)) is too small (must be > 3 cells on all axes) for a 6th order scheme to work...",
    )
  end

  cache = (
    xᵢ₊½=zeros(T, cellsize),
    ∂²x=zeros(T, cellsize),
    ∂x=zeros(T, cellsize),
    metric=zeros(T, cellsize),
    inner_deriv1=zeros(T, cellsize),
    outer_deriv1=zeros(T, cellsize),
    inner_deriv2=zeros(T, cellsize),
    outer_deriv2=zeros(T, cellsize),
  )

  return MonotoneExplicit6thOrderDiscretization(cache)
end

function ∂x∂ξ!(m::MonotoneExplicit6thOrderDiscretization, ∂x_∂ξ, x, domain, axis)
  xᵢ₊½ = m.cache.xᵢ₊½
  ∂²x = m.cache.∂²x
  ∂x = m.cache.∂x

  toedge!(xᵢ₊½, ∂²x, ∂x, x, domain, axis)

  inner_domain = expand(domain, axis, 0)
  for i in inner_domain
    ᵢ₋₁ = down(i, axis, 1)
    ∂x_∂ξ[i] = xᵢ₊½[i] - xᵢ₊½[ᵢ₋₁]
  end

  return nothing
end

function conserved_metric!(
  m::MonotoneExplicit6thOrderDiscretization, ξ̂x, y, η_axis, z, ζ_axis, domain
)

  # 1st term
  yη = m.cache.inner_deriv1
  yηz_ζ = m.cache.outer_deriv1

  ∂x∂ξ!(m, yη, y, domain, η_axis)
  yηz = mappedarray(*, yη, z)
  ∂x∂ξ!(m, yηz_ζ, yηz, domain, ζ_axis)

  # 2nd term
  yζ = m.cache.inner_deriv2
  yζz_η = m.cache.outer_deriv2

  ∂x∂ξ!(m, yζ, y, domain, ζ_axis) # -> yζ
  yζz = mappedarray(*, yζ, z)
  ∂x∂ξ!(m, yζz_η, yζz, domain, η_axis)

  for i in domain
    ξ̂x[i] = yηz_ζ[i] - yζz_η[i]
  end

  return nothing
end

function update_edge_conserved_metrics!(
  m::MonotoneExplicit6thOrderDiscretization,
  edge_metrics,
  centroid_coords::@NamedTuple{x::Array{T,3}, y::Array{T,3}, z::Array{T,3}},
  domain,
) where {T}
  x = centroid_coords.x
  y = centroid_coords.y
  z = centroid_coords.z

  # cache-arrays
  ∂²x = m.cache.∂²x
  ∂x = m.cache.∂x

  i₊½_domain = expand(domain, 1, 0)
  j₊½_domain = expand(domain, 2, 0)
  k₊½_domain = expand(domain, 3, 0)

  # ξ̂x = (y_η z)_ζ − (y_ζ z)_η
  ξ̂x = m.cache.metric
  conserved_metric!(m, ξ̂x, y, η, z, ζ, domain)
  toedge!(edge_metrics.i₊½.ξ̂x, ∂²x, ∂x, ξ̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂x, ∂²x, ∂x, ξ̂x, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ξ̂x, ∂²x, ∂x, ξ̂x, k₊½_domain, ζ)

  # η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
  η̂x = m.cache.metric
  conserved_metric!(m, η̂x, y, ζ, z, ξ, domain)
  toedge!(edge_metrics.i₊½.η̂x, ∂²x, ∂x, η̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂x, ∂²x, ∂x, η̂x, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.η̂x, ∂²x, ∂x, η̂x, k₊½_domain, ζ)

  # ζ̂x = (y_ξ z)_η − (y_η z)_ξ
  ζ̂x = m.cache.metric
  conserved_metric!(m, ζ̂x, y, ξ, z, η, domain)
  toedge!(edge_metrics.i₊½.ζ̂x, ∂²x, ∂x, ζ̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ζ̂x, ∂²x, ∂x, ζ̂x, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ζ̂x, ∂²x, ∂x, ζ̂x, k₊½_domain, ζ)

  # ξ̂y = (z_η x)_ζ − (z_ζ x)_η
  ξ̂y = m.cache.metric
  conserved_metric!(m, ξ̂y, z, η, x, ζ, domain)
  toedge!(edge_metrics.i₊½.ξ̂y, ∂²x, ∂x, ξ̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂y, ∂²x, ∂x, ξ̂y, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ξ̂y, ∂²x, ∂x, ξ̂y, k₊½_domain, ζ)

  # η̂y = (z_ζ x)_ξ − (z_ξ x)_ζ
  η̂y = m.cache.metric
  conserved_metric!(m, η̂y, z, ζ, x, ξ, domain)
  toedge!(edge_metrics.i₊½.η̂y, ∂²x, ∂x, η̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂y, ∂²x, ∂x, η̂y, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.η̂y, ∂²x, ∂x, η̂y, k₊½_domain, ζ)

  # ζ̂y = (z_ξ x)_η − (z_η x)_ξ
  ζ̂y = m.cache.metric
  conserved_metric!(m, ζ̂y, z, ξ, x, η, domain)
  toedge!(edge_metrics.i₊½.ζ̂y, ∂²x, ∂x, ζ̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ζ̂y, ∂²x, ∂x, ζ̂y, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ζ̂y, ∂²x, ∂x, ζ̂y, k₊½_domain, ζ)

  # ξ̂z = (x_η y)_ζ − (x_ζ y)_η
  ξ̂z = m.cache.metric
  conserved_metric!(m, ξ̂z, x, η, y, ζ, domain)
  toedge!(edge_metrics.i₊½.ξ̂z, ∂²x, ∂x, ξ̂z, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂z, ∂²x, ∂x, ξ̂z, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ξ̂z, ∂²x, ∂x, ξ̂z, k₊½_domain, ζ)

  # η̂z = (x_ζ y)_ξ − (x_ξ y)_ζ
  η̂z = m.cache.metric
  conserved_metric!(m, η̂z, x, ζ, y, ξ, domain)
  toedge!(edge_metrics.i₊½.η̂z, ∂²x, ∂x, η̂z, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂z, ∂²x, ∂x, η̂z, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.η̂z, ∂²x, ∂x, η̂z, k₊½_domain, ζ)

  # ζ̂z = (x_ξ y)_η − (x_η y)_ξ
  ζ̂z = m.cache.metric
  conserved_metric!(m, ζ̂z, x, ξ, y, η, domain)
  toedge!(edge_metrics.i₊½.ζ̂z, ∂²x, ∂x, ζ̂z, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ζ̂z, ∂²x, ∂x, ζ̂z, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ζ̂z, ∂²x, ∂x, ζ̂z, k₊½_domain, ζ)

  return nothing
end

function update_cell_center_metrics!(
  m::MonotoneExplicit6thOrderDiscretization,
  cell_center_metrics,
  centroid_coords::@NamedTuple{x::Array{T,3}, y::Array{T,3}, z::Array{T,3}},
  domain,
) where {T}
  xξ = @views cell_center_metrics.ξx
  yξ = @views cell_center_metrics.ξy
  zξ = @views cell_center_metrics.ξz
  xη = @views cell_center_metrics.ηx
  yη = @views cell_center_metrics.ηy
  zη = @views cell_center_metrics.ηz
  xζ = @views cell_center_metrics.ζx
  yζ = @views cell_center_metrics.ζy
  zζ = @views cell_center_metrics.ζz

  # First compute the inverse metrics (and store in the forward metric arrays for now)
  ∂x∂ξ!(m, cell_center_metrics.ξx, centroid_coords.x, domain, ξ) # xξ
  ∂x∂ξ!(m, cell_center_metrics.ξy, centroid_coords.y, domain, ξ) # yξ
  ∂x∂ξ!(m, cell_center_metrics.ξz, centroid_coords.z, domain, ξ) # zξ
  ∂x∂ξ!(m, cell_center_metrics.ηx, centroid_coords.x, domain, η) # xη
  ∂x∂ξ!(m, cell_center_metrics.ηy, centroid_coords.y, domain, η) # yη
  ∂x∂ξ!(m, cell_center_metrics.ηz, centroid_coords.z, domain, η) # zη
  ∂x∂ξ!(m, cell_center_metrics.ζx, centroid_coords.x, domain, ζ) # xζ
  ∂x∂ξ!(m, cell_center_metrics.ζy, centroid_coords.y, domain, ζ) # yζ
  ∂x∂ξ!(m, cell_center_metrics.ζz, centroid_coords.z, domain, ζ) # zζ

  # Now compute the jacobian matrix for each entry, and store the proper
  # metric back in it's place
  for idx in domain
    _jacobian = @SMatrix [
      xξ[idx] xη[idx] xζ[idx]
      yξ[idx] yη[idx] yζ[idx]
      zξ[idx] zη[idx] zζ[idx]
    ]

    cell_center_metrics.J[idx] = det(_jacobian)

    _inv_jacobian = inv(_jacobian)
    cell_center_metrics.ξx[idx] = _inv_jacobian[1, 1]
    cell_center_metrics.ξy[idx] = _inv_jacobian[1, 2]
    cell_center_metrics.ξz[idx] = _inv_jacobian[1, 3]
    cell_center_metrics.ηx[idx] = _inv_jacobian[2, 1]
    cell_center_metrics.ηy[idx] = _inv_jacobian[2, 2]
    cell_center_metrics.ηz[idx] = _inv_jacobian[2, 3]
    cell_center_metrics.ζx[idx] = _inv_jacobian[3, 1]
    cell_center_metrics.ζy[idx] = _inv_jacobian[3, 2]
    cell_center_metrics.ζz[idx] = _inv_jacobian[3, 3]
  end

  return nothing
end

end
