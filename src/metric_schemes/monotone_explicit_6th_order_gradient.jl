
include("operators.jl")

"""
A Monotone Explicit 6th Order Gradient (MEG6) scheme to discretize
the grid metrics
"""
struct MEG6Scheme{N,T,AA<:AbstractArray{T,N},B,C,D,E,F} <: AbstractMetricScheme
  ∂_cache::AA # cache array for derivative terms of any scalar term ϕ 
  ∂²_cache::AA # cache array for 2nd derivative terms of any scalar term ϕ 
  Jᵢ::AA # cell-centered Jacobian (det(jacobian matrix))
  ∂x∂ξᵢ::B # cell-centered inverse metrics ∂x/∂(ξ,η,ζ), ∂y/∂(ξ,η,ζ), ∂z/∂(ξ,η,ζ)
  ∂ξ∂xᵢ::C # cell-centered metrics ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  ∂ξ̂∂xᵢ::D # cell-centered conservative metrics, where ξ̂ = ξ/J; ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  ∂ξ̂∂xᵢ₊½::E # conservative metrics at the cell edges ((i,j,k)+1/2), where ξ̂ = ξ/J; ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  edge_cache::F
end

# MEG6Scheme Constructor
function MEG6Scheme(celldims::NTuple{N,Int}; T=Float64) where {N}
  ∂ϕᵢ = zeros(T, celldims)
  Jᵢ = zeros(T, celldims)
  ∂²ϕᵢ = zeros(T, celldims)

  ∂x∂ξᵢ = NamedTuple{(x_names(celldims))}( # (∂x, ∂y, ∂z)
    nested_tuple(celldims, ξ_names, T) for i in 1:N  # (∂ξ, ∂η, ∂ζ)
  )

  ∂ξ∂xᵢ = NamedTuple{(ξ_names(celldims))}( # (∂ξ, ∂η, ∂ζ)
    nested_tuple(celldims, x_names, T) for i in 1:N # (∂x, ∂y, ∂z)
  )

  ∂ξ̂∂xᵢ = NamedTuple{(ξ̂_names(celldims))}( # (∂ξ, ∂η, ∂ζ)
    nested_tuple(celldims, x_names, T) for i in 1:N # (∂x, ∂y, ∂z)
  )

  ∂ξ̂∂xᵢ₊½ = NamedTuple{(ξ_names(celldims))}( # (∂ξ, ∂η, ∂ζ)
    NamedTuple{(x_names(celldims))}( # (∂x, ∂y, ∂z)
      nested_tuple(celldims, edge_names, T) for i in 1:N # (i₊½, j₊½, k₊½)
    ) for j in 1:N
  )

  if N == 3
    # 3D needs another cache array for edge interpolation
    edge_cache = (zeros(T, celldims), zeros(T, celldims))
  else
    edge_cache = nothing
  end

  return MEG6Scheme(∂ϕᵢ, ∂²ϕᵢ, Jᵢ, ∂x∂ξᵢ, ∂ξ∂xᵢ, ∂ξ̂∂xᵢ, ∂ξ̂∂xᵢ₊½, edge_cache)
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

function update_metrics!(m::MEG6Scheme{3}, (xc, yc, zc), cell_center_metrics, edge_metrics)
  update_cell_center_metrics!(m, (xc, yc, zc), cell_center_metrics, edge_metrics)
  return update_edge_metrics!(m, (xc, yc, zc), cell_center_metrics, edge_metrics)
end

function update_cell_center_metrics!(
  ::MEG6Scheme{3}, (xc, yc, zc), cell_center_metrics, edge_metrics
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
  inner∂!(cell_center_metrics.ξx, x, ξ) # xξ
  inner∂!(cell_center_metrics.ξy, y, ξ) # yξ
  inner∂!(cell_center_metrics.ξz, z, ξ) # zξ
  inner∂!(cell_center_metrics.ηx, x, η) # xη
  inner∂!(cell_center_metrics.ηy, y, η) # yη
  inner∂!(cell_center_metrics.ηz, z, η) # zη
  inner∂!(cell_center_metrics.ζx, x, ζ) # xζ
  inner∂!(cell_center_metrics.ζy, y, ζ) # yζ
  inner∂!(cell_center_metrics.ζz, z, ζ) # zζ

  domain = CartesianIndices(x)

  # Now compute the jacobian matrix for each entry, and store the proper
  # metric back in it's place
  for idx in domain
    _jacobian = @SMatrix [
      xξ[idx] xη[idx] xζ[idx]
      yξ[idx] yη[idx] yζ[idx]
      zξ[idx] zη[idx] zζ[idx]
    ]

    _inv_jacobian = inv(_jacobian)
    J[idx] = det(_jacobian)

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

function update_edge_metrics!(
  m::MEG6Scheme{3}, (xc, yc, zc), cell_center_metrics, edge_metrics
)
  ξ, η, ζ = (1, 2, 3)

  to_edge!(edge_metrics.ᵢ₊½.J, cell_center_metrics.J, ξ)
  to_edge!(edge_metrics.ⱼ₊½.J, cell_center_metrics.J, η)
  to_edge!(edge_metrics.ₖ₊½.J, cell_center_metrics.J, ζ)

  ξ̂_x = m.ξ̂cache
  conserved_metric(ξ̂_x, (yc, zc), (η, ζ), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.ξ̂x, ξ̂_x, ξ)
  to_edge!(edge_metrics.ⱼ₊½.ξ̂x, ξ̂_x, η)
  to_edge!(edge_metrics.ₖ₊½.ξ̂x, ξ̂_x, ζ)

  η̂_x = m.ξ̂cache
  conserved_metric(η̂_x, (yc, zc), (ζ, ξ), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.η̂x, η̂_x, ξ)
  to_edge!(edge_metrics.ⱼ₊½.η̂x, η̂_x, η)
  to_edge!(edge_metrics.ₖ₊½.η̂x, η̂_x, ζ)

  ζ̂_x = m.ξ̂cache
  conserved_metric(ζ̂_x, (yc, zc), (ξ, η), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.ζ̂x, ζ̂_x, ξ)
  to_edge!(edge_metrics.ⱼ₊½.ζ̂x, ζ̂_x, η)
  to_edge!(edge_metrics.ₖ₊½.ζ̂x, ζ̂_x, ζ)

  ξ̂_y = m.ξ̂cache
  conserved_metric(ξ̂_y, (zc, xc), (η, ζ), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.ξ̂y, ξ̂_y, ξ)
  to_edge!(edge_metrics.ⱼ₊½.ξ̂y, ξ̂_y, η)
  to_edge!(edge_metrics.ₖ₊½.ξ̂y, ξ̂_y, ζ)

  η̂_y = m.ξ̂cache
  conserved_metric(η̂_y, (zc, xc), (ζ, ξ), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.η̂y, η̂_y, ξ)
  to_edge!(edge_metrics.ⱼ₊½.η̂y, η̂_y, η)
  to_edge!(edge_metrics.ₖ₊½.η̂y, η̂_y, ζ)

  ζ̂_y = m.ξ̂cache
  conserved_metric(ζ̂_y, (zc, xc), (ξ, η), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.ζ̂y, ζ̂_y, ξ)
  to_edge!(edge_metrics.ⱼ₊½.ζ̂y, ζ̂_y, η)
  to_edge!(edge_metrics.ₖ₊½.ζ̂y, ζ̂_y, ζ)

  ξ̂_z = m.ξ̂cache
  conserved_metric(ξ̂_z, (xc, yc), (η, ζ), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.ξ̂z, ξ̂_z, ξ)
  to_edge!(edge_metrics.ⱼ₊½.ξ̂z, ξ̂_z, η)
  to_edge!(edge_metrics.ₖ₊½.ξ̂z, ξ̂_z, ζ)

  η̂_z = m.ξ̂cache
  conserved_metric(η̂_z, (xc, yc), (ζ, ξ), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.η̂z, η̂_z, ξ)
  to_edge!(edge_metrics.ⱼ₊½.η̂z, η̂_z, η)
  to_edge!(edge_metrics.ₖ₊½.η̂z, η̂_z, ζ)

  ζ̂_z = m.ξ̂cache
  conserved_metric(ζ̂_z, (xc, yc), (ξ, η), m.cache)
  to_edge!(edge_metrics.ᵢ₊½.ζ̂z, ζ̂_z, ξ)
  to_edge!(edge_metrics.ⱼ₊½.ζ̂z, ζ̂_z, η)
  to_edge!(edge_metrics.ₖ₊½.ζ̂z, ζ̂_z, ζ)

  return nothing
end