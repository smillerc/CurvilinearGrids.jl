
include("operators.jl")

const ξ = 1
const η = 2
const ζ = 3

"""
A Monotone Explicit 6th Order Gradient (MEG6) scheme to discretize
the grid metrics
"""
struct MEG6Scheme{C} <: AbstractMetricScheme
  cache::C
end

# MEG6Scheme Constructor
function MEG6Scheme(celldims::NTuple{N,Int}; T=Float64) where {N}
  cache = (
    ξ̂=zeros(T, celldims), mixed_deriv1=zeros(T, celldims), mixed_deriv2=zeros(T, celldims)
  )
  return MEG6Scheme(cache)
end

# Inner partial product, i.e. (∂x/∂ξ)z
function inner∂prod!(
  yξ_z::AbstractArray{T,N}, y::AbstractArray{T,N}, z::AbstractArray{T,N}, ξ::Int
) where {T,N}
  f1 = 2282 / 2880
  f2 = -546 / 2880
  f3 = 89 / 2880
  f4 = -3 / 2880
  f5 = -1 / 2880

  domain = CartesianIndices(y)
  inner_domain = expand(domain, -5)

  @batch for i in inner_domain
    ᵢ₊₁ = up(i, ξ, 1)
    ᵢ₊₂ = up(i, ξ, 2)
    ᵢ₊₃ = up(i, ξ, 3)
    ᵢ₊₄ = up(i, ξ, 4)
    ᵢ₊₅ = up(i, ξ, 5)
    ᵢ₋₁ = down(i, ξ, 1)
    ᵢ₋₂ = down(i, ξ, 2)
    ᵢ₋₃ = down(i, ξ, 3)
    ᵢ₋₄ = down(i, ξ, 4)
    ᵢ₋₅ = down(i, ξ, 5)

    yξ_z[i] =
      (
        f1 * (y[ᵢ₊₁] - y[ᵢ₋₁]) + #
        f2 * (y[ᵢ₊₂] - y[ᵢ₋₂]) + #
        f3 * (y[ᵢ₊₃] - y[ᵢ₋₃]) + #
        f4 * (y[ᵢ₊₄] - y[ᵢ₋₄]) + #
        f5 * (y[ᵢ₊₅] - y[ᵢ₋₅])
      ) * z[i]
  end
end

# Inner partial, i.e. (∂x/∂ξ)
# TODO: make the boundary versions of this - expand the 
function inner∂!(yξ::AbstractArray{T,N}, y::AbstractArray{T,N}, ξ::Int) where {T,N}
  f1 = 2282 / 2880
  f2 = -546 / 2880
  f3 = 89 / 2880
  f4 = -3 / 2880
  f5 = -1 / 2880

  domain = CartesianIndices(y)
  inner_domain = expand(domain, -5)

  @batch for i in inner_domain
    ᵢ₊₁ = up(i, ξ, 1)
    ᵢ₊₂ = up(i, ξ, 2)
    ᵢ₊₃ = up(i, ξ, 3)
    ᵢ₊₄ = up(i, ξ, 4)
    ᵢ₊₅ = up(i, ξ, 5)
    ᵢ₋₁ = down(i, ξ, 1)
    ᵢ₋₂ = down(i, ξ, 2)
    ᵢ₋₃ = down(i, ξ, 3)
    ᵢ₋₄ = down(i, ξ, 4)
    ᵢ₋₅ = down(i, ξ, 5)

    yξ[i] =
      f1 * (y[ᵢ₊₁] - y[ᵢ₋₁]) + #
      f2 * (y[ᵢ₊₂] - y[ᵢ₋₂]) + #
      f3 * (y[ᵢ₊₃] - y[ᵢ₋₃]) + #
      f4 * (y[ᵢ₊₄] - y[ᵢ₋₄]) + #
      f5 * (y[ᵢ₊₅] - y[ᵢ₋₅])
  end
end

function conserved!(
  ξ̂x::AbstractArray{T,N}, # conserved metric, could be ξ̂x, η̂z, ζ̂y, etc..
  y_ηz::AbstractArray{T,N}, # 1st term, (∂y/∂η)⋅z, in the 2-term eq to determine ξ̂x
  ζ::Int, # outer dimension to differentiate the 1st term on, e.g. ξ:i±1/2, η:j±1/2, ζ:k±1/2
  yζ_z::AbstractArray{T,N}, # 2nd term, (∂y/∂ζ)⋅z, in the 2-term eq to determine ξ̂x
  η::Int, # outer dimension to differentiate the 2nd term on, e.g. ξ:i±1/2, η:j±1/2, ζ:k±1/2
) where {T,N}
  a = 2282 / 2880
  b = -546 / 2880
  c = 89 / 2880
  d = -3 / 2880
  e = -1 / 2880

  domain = CartesianIndices(ξ̂x)
  inner_domain = expand(domain, -5) # shrink due to the stencil size

  @batch for i in inner_domain
    ₖ₊₁ = up(i, ζ, 1)
    ₖ₊₂ = up(i, ζ, 2)
    ₖ₊₃ = up(i, ζ, 3)
    ₖ₊₄ = up(i, ζ, 4)
    ₖ₊₅ = up(i, ζ, 5)
    ₖ₋₁ = down(i, ζ, 1)
    ₖ₋₂ = down(i, ζ, 2)
    ₖ₋₃ = down(i, ζ, 3)
    ₖ₋₄ = down(i, ζ, 4)
    ₖ₋₅ = down(i, ζ, 5)

    ⱼ₊₁ = up(i, η, 1)
    ⱼ₊₂ = up(i, η, 2)
    ⱼ₊₃ = up(i, η, 3)
    ⱼ₊₄ = up(i, η, 4)
    ⱼ₊₅ = up(i, η, 5)
    ⱼ₋₁ = down(i, η, 1)
    ⱼ₋₂ = down(i, η, 2)
    ⱼ₋₃ = down(i, η, 3)
    ⱼ₋₄ = down(i, η, 4)
    ⱼ₋₅ = down(i, η, 5)

    ξ̂x[i] =
      (
        a * (y_ηz[ₖ₊₁] - y_ηz[ₖ₋₁]) +
        b * (y_ηz[ₖ₊₂] - y_ηz[ₖ₋₂]) +
        c * (y_ηz[ₖ₊₃] - y_ηz[ₖ₋₃]) +
        d * (y_ηz[ₖ₊₄] - y_ηz[ₖ₋₄]) +
        e * (y_ηz[ₖ₊₅] - y_ηz[ₖ₋₅])
      ) - (
        a * (yζ_z[ⱼ₊₁] - yζ_z[ⱼ₋₁]) +
        b * (yζ_z[ⱼ₊₂] - yζ_z[ⱼ₋₂]) +
        c * (yζ_z[ⱼ₊₃] - yζ_z[ⱼ₋₃]) +
        d * (yζ_z[ⱼ₊₄] - yζ_z[ⱼ₋₄]) +
        e * (yζ_z[ⱼ₊₅] - yζ_z[ⱼ₋₅])
      )
  end

  return nothing
end

function to_edge!(ϕᵢ₊½::AbstractArray{T,N}, ϕ::AbstractArray{T,N}, dim::Int) where {T,N}
  a = 1821 / 2880
  b = -461 / 2880
  c = 85 / 2880
  d = -4 / 2880
  e = -1 / 2880

  domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, -5) # shrink due to the stencil size

  @batch for i in inner_domain
    ᵢ₊₁ = up(i, dim, 1)
    ᵢ₊₂ = up(i, dim, 2)
    ᵢ₊₃ = up(i, dim, 3)
    ᵢ₊₄ = up(i, dim, 4)
    ᵢ₊₅ = up(i, dim, 5)
    ᵢ₋₁ = down(i, dim, 1)
    ᵢ₋₂ = down(i, dim, 2)
    ᵢ₋₃ = down(i, dim, 3)
    ᵢ₋₄ = down(i, dim, 4)

    ϕᵢ₊½[i] = (
      a * (ϕ[ᵢ₊₁] + ϕ[i]) +
      b * (ϕ[ᵢ₊₂] + ϕ[ᵢ₋₁]) +
      c * (ϕ[ᵢ₊₃] + ϕ[ᵢ₋₂]) +
      d * (ϕ[ᵢ₊₄] + ϕ[ᵢ₋₃]) +
      e * (ϕ[ᵢ₊₅] + ϕ[ᵢ₋₄])
    )
  end

  return nothing
end

function conserved_metric(
  ξ̂_x::AbstractArray{T,N}, (y, z), (η, ζ)::NTuple{2,Int}, cache
) where {T,N}
  yη_z = @views cache.mixed_deriv1
  yζ_z = @views cache.mixed_deriv2

  inner∂prod!(yη_z, y, z, η) # -> yη_z
  inner∂prod!(yζ_z, y, z, ζ) # -> yζ_z
  conserved!(ξ̂_x, yη_z, ζ, yζ_z, η) # -> ξ̂x

  return ξ̂_x
end

function update_metrics_3d!(
  m::MEG6Scheme, centroid_coords, cell_center_metrics, edge_metrics
)
  update_cell_center_metrics_3d!(m, centroid_coords, cell_center_metrics, edge_metrics)
  update_edge_metrics_3d!(m, centroid_coords, cell_center_metrics, edge_metrics)
  return nothing
end

function update_cell_center_metrics_3d!(
  ::MEG6Scheme, centroid_coords, cell_center_metrics, edge_metrics
)
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
  inner∂!(cell_center_metrics.ξx, centroid_coords.x, ξ) # xξ
  inner∂!(cell_center_metrics.ξy, centroid_coords.y, ξ) # yξ
  inner∂!(cell_center_metrics.ξz, centroid_coords.z, ξ) # zξ
  inner∂!(cell_center_metrics.ηx, centroid_coords.x, η) # xη
  inner∂!(cell_center_metrics.ηy, centroid_coords.y, η) # yη
  inner∂!(cell_center_metrics.ηz, centroid_coords.z, η) # zη
  inner∂!(cell_center_metrics.ζx, centroid_coords.x, ζ) # xζ
  inner∂!(cell_center_metrics.ζy, centroid_coords.y, ζ) # yζ
  inner∂!(cell_center_metrics.ζz, centroid_coords.z, ζ) # zζ

  domain = CartesianIndices(centroid_coords.x)

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

function update_edge_metrics_3d!(
  m::MEG6Scheme, centroid_coords, cell_center_metrics, edge_metrics
)
  xc = centroid_coords.x
  yc = centroid_coords.y
  zc = centroid_coords.z

  to_edge!(edge_metrics.i₊½.J, cell_center_metrics.J, ξ)
  to_edge!(edge_metrics.j₊½.J, cell_center_metrics.J, η)
  to_edge!(edge_metrics.k₊½.J, cell_center_metrics.J, ζ)

  ξ̂_x = m.cache.ξ̂
  conserved_metric(ξ̂_x, (yc, zc), (η, ζ), m.cache)
  to_edge!(edge_metrics.i₊½.ξ̂x, ξ̂_x, ξ)
  to_edge!(edge_metrics.j₊½.ξ̂x, ξ̂_x, η)
  to_edge!(edge_metrics.k₊½.ξ̂x, ξ̂_x, ζ)

  η̂_x = m.cache.ξ̂
  conserved_metric(η̂_x, (yc, zc), (ζ, ξ), m.cache)
  to_edge!(edge_metrics.i₊½.η̂x, η̂_x, ξ)
  to_edge!(edge_metrics.j₊½.η̂x, η̂_x, η)
  to_edge!(edge_metrics.k₊½.η̂x, η̂_x, ζ)

  ζ̂_x = m.cache.ξ̂
  conserved_metric(ζ̂_x, (yc, zc), (ξ, η), m.cache)
  to_edge!(edge_metrics.i₊½.ζ̂x, ζ̂_x, ξ)
  to_edge!(edge_metrics.j₊½.ζ̂x, ζ̂_x, η)
  to_edge!(edge_metrics.k₊½.ζ̂x, ζ̂_x, ζ)

  ξ̂_y = m.cache.ξ̂
  conserved_metric(ξ̂_y, (zc, xc), (η, ζ), m.cache)
  to_edge!(edge_metrics.i₊½.ξ̂y, ξ̂_y, ξ)
  to_edge!(edge_metrics.j₊½.ξ̂y, ξ̂_y, η)
  to_edge!(edge_metrics.k₊½.ξ̂y, ξ̂_y, ζ)

  η̂_y = m.cache.ξ̂
  conserved_metric(η̂_y, (zc, xc), (ζ, ξ), m.cache)
  to_edge!(edge_metrics.i₊½.η̂y, η̂_y, ξ)
  to_edge!(edge_metrics.j₊½.η̂y, η̂_y, η)
  to_edge!(edge_metrics.k₊½.η̂y, η̂_y, ζ)

  ζ̂_y = m.cache.ξ̂
  conserved_metric(ζ̂_y, (zc, xc), (ξ, η), m.cache)
  to_edge!(edge_metrics.i₊½.ζ̂y, ζ̂_y, ξ)
  to_edge!(edge_metrics.j₊½.ζ̂y, ζ̂_y, η)
  to_edge!(edge_metrics.k₊½.ζ̂y, ζ̂_y, ζ)

  ξ̂_z = m.cache.ξ̂
  conserved_metric(ξ̂_z, (xc, yc), (η, ζ), m.cache)
  to_edge!(edge_metrics.i₊½.ξ̂z, ξ̂_z, ξ)
  to_edge!(edge_metrics.j₊½.ξ̂z, ξ̂_z, η)
  to_edge!(edge_metrics.k₊½.ξ̂z, ξ̂_z, ζ)

  η̂_z = m.cache.ξ̂
  conserved_metric(η̂_z, (xc, yc), (ζ, ξ), m.cache)
  to_edge!(edge_metrics.i₊½.η̂z, η̂_z, ξ)
  to_edge!(edge_metrics.j₊½.η̂z, η̂_z, η)
  to_edge!(edge_metrics.k₊½.η̂z, η̂_z, ζ)

  ζ̂_z = m.cache.ξ̂
  conserved_metric(ζ̂_z, (xc, yc), (ξ, η), m.cache)
  to_edge!(edge_metrics.i₊½.ζ̂z, ζ̂_z, ξ)
  to_edge!(edge_metrics.j₊½.ζ̂z, ζ̂_z, η)
  to_edge!(edge_metrics.k₊½.ζ̂z, ζ̂_z, ζ)

  return nothing
end