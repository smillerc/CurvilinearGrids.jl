module MonotoneExplicit6thOrderScheme

using LinearAlgebra
using StaticArrays, MappedArrays, StructArrays
using Polyester
using KernelAbstractions

using ..IndexingUtils

export MonotoneExplicit6thOrderDiscretization
export update_edge_conserved_metrics!, update_cell_center_metrics!

const ξ = 1
const η = 2
const ζ = 3

include("gradients.jl")
include("interpolation.jl")

CacheType{T,N,AA<:AbstractArray{T,N}} = @NamedTuple{
  xᵢ₊½::AA,
  ∂²x::AA,
  ∂x::AA,
  metric::AA,
  inner_deriv1::AA,
  outer_deriv1::AA,
  inner_deriv2::AA,
  outer_deriv2::AA,
}

# struct MonotoneExplicit6thOrderDiscretization{T,N,AA}
#   cache::CacheType{T,N,AA}
# end
struct MonotoneExplicit6thOrderDiscretization{T,N,AA<:AbstractArray{T,N}}
  cache::@NamedTuple{
    xᵢ₊½::AA,
    ∂²x::AA,
    ∂x::AA,
    metric::AA,
    inner_deriv1::AA,
    outer_deriv1::AA,
    inner_deriv2::AA,
    outer_deriv2::AA,
  }
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

  AA = typeof(cache.metric)
  return MonotoneExplicit6thOrderDiscretization{T,N,AA}(cache)
end

function ∂x∂ξ!(
  m::MonotoneExplicit6thOrderDiscretization,
  ∂x_∂ξ::AbstractArray{T},
  x,
  domain,
  axis,
  ϵ=10eps(T),
) where {T}
  xᵢ₊½ = m.cache.xᵢ₊½
  ∂²x = m.cache.∂²x
  ∂x = m.cache.∂x

  toedge!(xᵢ₊½, ∂²x, ∂x, x, domain, axis)

  inner_domain = expand(domain, axis, 0)
  @batch for i in inner_domain
    ᵢ₋₁ = down(i, axis, 1)
    # ∂x_∂ξ[i] = xᵢ₊½[i] - xᵢ₊½[ᵢ₋₁]
    # ∂x_∂ξ[i] = xᵢ₊½[i] - xᵢ₊½[ᵢ₋₁]
    _∂x = xᵢ₊½[i] - xᵢ₊½[ᵢ₋₁]

    ∂x_∂ξ[i] = _∂x * (abs(_∂x) >= ϵ)
  end

  return nothing
end

function conserved_metric!(
  m::MonotoneExplicit6thOrderDiscretization,
  ξ̂x::AbstractArray{T},
  y,
  η_axis,
  z,
  ζ_axis,
  domain,
  ϵ=10eps(T),
) where {T}

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

  @batch for i in domain
    # ξ̂x[i] = yηz_ζ[i] - yζz_η[i]
    _ξ̂x = yηz_ζ[i] - yζz_η[i]
    ξ̂x[i] = _ξ̂x * (abs(_ξ̂x) >= ϵ)
  end

  return nothing
end

function update_edge_conserved_metrics!(
  m::MonotoneExplicit6thOrderDiscretization{T,3},
  edge_metrics,
  cell_center_metrics,
  centroid_coords,
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
  toedge!(edge_metrics.i₊½.ξ̂.x₁, ∂²x, ∂x, ξ̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂.x₁, ∂²x, ∂x, ξ̂x, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ξ̂.x₁, ∂²x, ∂x, ξ̂x, k₊½_domain, ζ)

  # η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
  η̂x = m.cache.metric
  conserved_metric!(m, η̂x, y, ζ, z, ξ, domain)
  toedge!(edge_metrics.i₊½.η̂.x₁, ∂²x, ∂x, η̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂.x₁, ∂²x, ∂x, η̂x, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.η̂.x₁, ∂²x, ∂x, η̂x, k₊½_domain, ζ)

  # ζ̂x = (y_ξ z)_η − (y_η z)_ξ
  ζ̂x = m.cache.metric
  conserved_metric!(m, ζ̂x, y, ξ, z, η, domain)
  toedge!(edge_metrics.i₊½.ζ̂.x₁, ∂²x, ∂x, ζ̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ζ̂.x₁, ∂²x, ∂x, ζ̂x, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ζ̂.x₁, ∂²x, ∂x, ζ̂x, k₊½_domain, ζ)

  # ξ̂y = (z_η x)_ζ − (z_ζ x)_η
  ξ̂y = m.cache.metric
  conserved_metric!(m, ξ̂y, z, η, x, ζ, domain)
  toedge!(edge_metrics.i₊½.ξ̂.x₂, ∂²x, ∂x, ξ̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂.x₂, ∂²x, ∂x, ξ̂y, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ξ̂.x₂, ∂²x, ∂x, ξ̂y, k₊½_domain, ζ)

  # η̂y = (z_ζ x)_ξ − (z_ξ x)_ζ
  η̂y = m.cache.metric
  conserved_metric!(m, η̂y, z, ζ, x, ξ, domain)
  toedge!(edge_metrics.i₊½.η̂.x₂, ∂²x, ∂x, η̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂.x₂, ∂²x, ∂x, η̂y, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.η̂.x₂, ∂²x, ∂x, η̂y, k₊½_domain, ζ)

  # ζ̂y = (z_ξ x)_η − (z_η x)_ξ
  ζ̂y = m.cache.metric
  conserved_metric!(m, ζ̂y, z, ξ, x, η, domain)
  toedge!(edge_metrics.i₊½.ζ̂.x₂, ∂²x, ∂x, ζ̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ζ̂.x₂, ∂²x, ∂x, ζ̂y, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ζ̂.x₂, ∂²x, ∂x, ζ̂y, k₊½_domain, ζ)

  # ξ̂z = (x_η y)_ζ − (x_ζ y)_η
  ξ̂z = m.cache.metric
  conserved_metric!(m, ξ̂z, x, η, y, ζ, domain)
  toedge!(edge_metrics.i₊½.ξ̂.x₃, ∂²x, ∂x, ξ̂z, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂.x₃, ∂²x, ∂x, ξ̂z, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ξ̂.x₃, ∂²x, ∂x, ξ̂z, k₊½_domain, ζ)

  # η̂z = (x_ζ y)_ξ − (x_ξ y)_ζ
  η̂z = m.cache.metric
  conserved_metric!(m, η̂z, x, ζ, y, ξ, domain)
  toedge!(edge_metrics.i₊½.η̂.x₃, ∂²x, ∂x, η̂z, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂.x₃, ∂²x, ∂x, η̂z, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.η̂.x₃, ∂²x, ∂x, η̂z, k₊½_domain, ζ)

  # ζ̂z = (x_ξ y)_η − (x_η y)_ξ
  ζ̂z = m.cache.metric
  conserved_metric!(m, ζ̂z, x, ξ, y, η, domain)
  toedge!(edge_metrics.i₊½.ζ̂.x₃, ∂²x, ∂x, ζ̂z, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ζ̂.x₃, ∂²x, ∂x, ζ̂z, j₊½_domain, η)
  toedge!(edge_metrics.k₊½.ζ̂.x₃, ∂²x, ∂x, ζ̂z, k₊½_domain, ζ)

  return nothing
end

function update_edge_conserved_metrics!(
  m::MonotoneExplicit6thOrderDiscretization{T,2},
  edge_metrics,
  cell_center_metrics,
  centroid_coords,
  domain,
) where {T}

  # cache-arrays
  ∂²x = m.cache.∂²x
  ∂x = m.cache.∂x

  i₊½_domain = expand(domain, 1, 0)
  j₊½_domain = expand(domain, 2, 0)

  ξ̂x = mappedarray(*, cell_center_metrics.ξ.x₁, cell_center_metrics.J)
  toedge!(edge_metrics.i₊½.ξ̂.x₁, ∂²x, ∂x, ξ̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂.x₁, ∂²x, ∂x, ξ̂x, j₊½_domain, η)

  # η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
  η̂x = mappedarray(*, cell_center_metrics.η.x₁, cell_center_metrics.J)
  toedge!(edge_metrics.i₊½.η̂.x₁, ∂²x, ∂x, η̂x, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂.x₁, ∂²x, ∂x, η̂x, j₊½_domain, η)

  ξ̂y = mappedarray(*, cell_center_metrics.ξ.x₂, cell_center_metrics.J)
  toedge!(edge_metrics.i₊½.ξ̂.x₂, ∂²x, ∂x, ξ̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.ξ̂.x₂, ∂²x, ∂x, ξ̂y, j₊½_domain, η)

  η̂y = mappedarray(*, cell_center_metrics.η.x₂, cell_center_metrics.J)
  toedge!(edge_metrics.i₊½.η̂.x₂, ∂²x, ∂x, η̂y, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.η̂.x₂, ∂²x, ∂x, η̂y, j₊½_domain, η)

  toedge!(edge_metrics.i₊½.J, ∂²x, ∂x, cell_center_metrics.J, i₊½_domain, ξ)
  toedge!(edge_metrics.j₊½.J, ∂²x, ∂x, cell_center_metrics.J, j₊½_domain, η)

  return nothing
end

function update_edge_conserved_metrics!(
  m::MonotoneExplicit6thOrderDiscretization{T,1},
  edge_metrics,
  cell_center_metrics,
  centroid_coords,
  domain,
) where {T}

  # cache-arrays
  ∂²x = m.cache.∂²x
  ∂x = m.cache.∂x

  i₊½_domain = expand(domain, 1, 0)

  ξ̂x = mappedarray(*, cell_center_metrics.ξ.x₁, cell_center_metrics.J)

  toedge!(edge_metrics.i₊½.ξ̂.x₁, ∂²x, ∂x, ξ̂x, i₊½_domain, ξ)

  toedge!(edge_metrics.i₊½.J, ∂²x, ∂x, cell_center_metrics.J, i₊½_domain, ξ)

  return nothing
end

function update_cell_center_metrics!(
  m::MonotoneExplicit6thOrderDiscretization{T,3},
  cell_center_metrics,
  centroid_coords,
  domain,
) where {T}
  xξ = cell_center_metrics.x₁.ξ
  yξ = cell_center_metrics.x₂.ξ
  zξ = cell_center_metrics.x₃.ξ
  xη = cell_center_metrics.x₁.η
  yη = cell_center_metrics.x₂.η
  zη = cell_center_metrics.x₃.η
  xζ = cell_center_metrics.x₁.ζ
  yζ = cell_center_metrics.x₂.ζ
  zζ = cell_center_metrics.x₃.ζ

  # First compute the inverse metrics (and store in the forward metric arrays for now)
  ∂x∂ξ!(m, xξ, centroid_coords.x, domain, ξ) # xξ
  ∂x∂ξ!(m, yξ, centroid_coords.y, domain, ξ) # yξ
  ∂x∂ξ!(m, zξ, centroid_coords.z, domain, ξ) # zξ
  ∂x∂ξ!(m, xη, centroid_coords.x, domain, η) # xη
  ∂x∂ξ!(m, yη, centroid_coords.y, domain, η) # yη
  ∂x∂ξ!(m, zη, centroid_coords.z, domain, η) # zη
  ∂x∂ξ!(m, xζ, centroid_coords.x, domain, ζ) # xζ
  ∂x∂ξ!(m, yζ, centroid_coords.y, domain, ζ) # yζ
  ∂x∂ξ!(m, zζ, centroid_coords.z, domain, ζ) # zζ

  # Now compute the jacobian matrix for each entry, and store the proper
  # metric back in it's place
  @batch for idx in domain
    _jacobian = @SMatrix [
      xξ[idx] xη[idx] xζ[idx]
      yξ[idx] yη[idx] yζ[idx]
      zξ[idx] zη[idx] zζ[idx]
    ]

    cell_center_metrics.J[idx] = det(_jacobian)

    _inv_jacobian = inv(_jacobian)
    cell_center_metrics.ξ.x₁[idx] = _inv_jacobian[1, 1]
    cell_center_metrics.ξ.x₂[idx] = _inv_jacobian[1, 2]
    cell_center_metrics.ξ.x₃[idx] = _inv_jacobian[1, 3]

    cell_center_metrics.η.x₁[idx] = _inv_jacobian[2, 1]
    cell_center_metrics.η.x₂[idx] = _inv_jacobian[2, 2]
    cell_center_metrics.η.x₃[idx] = _inv_jacobian[2, 3]

    cell_center_metrics.ζ.x₁[idx] = _inv_jacobian[3, 1]
    cell_center_metrics.ζ.x₂[idx] = _inv_jacobian[3, 2]
    cell_center_metrics.ζ.x₃[idx] = _inv_jacobian[3, 3]
  end

  return nothing
end

function update_cell_center_metrics!(
  m::MonotoneExplicit6thOrderDiscretization{T,2},
  cell_center_metrics,
  centroid_coords,
  domain,
) where {T}
  xξ = cell_center_metrics.x₁.ξ
  yξ = cell_center_metrics.x₂.ξ
  xη = cell_center_metrics.x₁.η
  yη = cell_center_metrics.x₂.η

  # First compute the inverse metrics (and store in the forward metric arrays for now)
  ∂x∂ξ!(m, xξ, centroid_coords.x, domain, ξ)
  ∂x∂ξ!(m, yξ, centroid_coords.y, domain, ξ)
  ∂x∂ξ!(m, xη, centroid_coords.x, domain, η)
  ∂x∂ξ!(m, yη, centroid_coords.y, domain, η)

  # Now compute the jacobian matrix for each entry, and store the proper
  # metric back in it's place
  @batch for idx in domain
    _jacobian = @SMatrix [
      xξ[idx] xη[idx]
      yξ[idx] yη[idx]
    ]

    cell_center_metrics.J[idx] = det(_jacobian)

    _inv_jacobian = inv(_jacobian)
    cell_center_metrics.ξ.x₁[idx] = _inv_jacobian[1, 1]
    cell_center_metrics.ξ.x₂[idx] = _inv_jacobian[1, 2]
    cell_center_metrics.η.x₁[idx] = _inv_jacobian[2, 1]
    cell_center_metrics.η.x₂[idx] = _inv_jacobian[2, 2]
  end

  return nothing
end

function update_cell_center_metrics!(
  m::MonotoneExplicit6thOrderDiscretization{T,1},
  cell_center_metrics,
  centroid_coords,
  domain,
) where {T}
  xξ = cell_center_metrics.x₁.ξ

  # First compute the inverse metrics (and store in the forward metric arrays for now)
  ∂x∂ξ!(m, xξ, centroid_coords.x, domain, ξ)

  # Now compute the jacobian matrix for each entry, and store the proper
  # metric back in it's place
  @batch for idx in domain
    _jacobian = xξ[idx]
    cell_center_metrics.J[idx] = _jacobian
    _inv_jacobian = inv(_jacobian)
    cell_center_metrics.ξ.x₁[idx] = _inv_jacobian
  end

  return nothing
end

end
