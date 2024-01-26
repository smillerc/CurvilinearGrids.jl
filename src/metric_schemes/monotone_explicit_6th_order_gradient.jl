
"""
A Monotone Explicit 6th Order Gradient (MEG6) scheme to discretize
the grid metrics
"""
struct MEG6Scheme{N,T,AA<:AbstractArray{T,N},B,C,D,E} <: AbstractMetricScheme
  ∂_cache::AA # cache array for derivative terms of any scalar term ϕ 
  ∂²_cache::AA # cache array for 2nd derivative terms of any scalar term ϕ 
  Jᵢ::AA # cell-centered Jacobian (det(jacobian matrix))
  ∂x∂ξᵢ::B # cell-centered inverse metrics ∂x/∂(ξ,η,ζ), ∂y/∂(ξ,η,ζ), ∂z/∂(ξ,η,ζ)
  ∂ξ∂xᵢ::C # cell-centered metrics ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  ∂ξ̂∂xᵢ::D # cell-centered conservative metrics, where ξ̂ = ξ/J; ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  ∂ξ̂∂xᵢ₊½::E # conservative metrics at the cell edges ((i,j,k)+1/2), where ξ̂ = ξ/J; ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
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
  return MEG6Scheme(∂ϕᵢ, ∂²ϕᵢ, Jᵢ, ∂x∂ξᵢ, ∂ξ∂xᵢ, ∂ξ̂∂xᵢ, ∂ξ̂∂xᵢ₊½)
end

# 3D version
"""
    update_metrics!(scheme, xc, yc, zc, domain) 

Update the 3D grid metrics given the cell-centroid positions `xc`, `yc`, and `zc`
"""
function update_metrics!(
  scheme::MEG6Scheme{3},
  xc::AbstractArray{T,3},
  yc::AbstractArray{T,3},
  zc::AbstractArray{T,3},
) where {T}
  # compute the ∂x/∂ξ|ᵢ, ∂x/∂η|ᵢ terms

  # Each dimension is on a different axes in the various metric arrays.
  # These are used since the metric functions are defined for
  # arbitrary dimensions, so we need to define what axis we're working on
  ξ_ax, η_ax, ζ_ax = (1, 2, 3)

  _cell_center_metric(scheme, xc, scheme.∂xᵢ∂ξᵢ.x.ξ, ξ_ax, domain) # ∂x/∂ξ
  _cell_center_metric(scheme, xc, scheme.∂xᵢ∂ξᵢ.x.η, η_ax, domain) # ∂x/∂η
  _cell_center_metric(scheme, xc, scheme.∂xᵢ∂ξᵢ.x.ζ, ζ_ax, domain) # ∂x/∂ζ

  _cell_center_metric(scheme, yc, scheme.∂xᵢ∂ξᵢ.y.ξ, ξ_ax, domain) # ∂y/∂ξ
  _cell_center_metric(scheme, yc, scheme.∂xᵢ∂ξᵢ.y.η, η_ax, domain) # ∂y/∂η
  _cell_center_metric(scheme, yc, scheme.∂xᵢ∂ξᵢ.y.ζ, ζ_ax, domain) # ∂y/∂ζ

  _cell_center_metric(scheme, zc, scheme.∂xᵢ∂ξᵢ.z.ξ, ξ_ax, domain) # ∂z/∂ξ
  _cell_center_metric(scheme, zc, scheme.∂xᵢ∂ξᵢ.z.η, η_ax, domain) # ∂z/∂η
  _cell_center_metric(scheme, zc, scheme.∂xᵢ∂ξᵢ.z.ζ, ζ_ax, domain) # ∂z/∂ζ

  # compute the inverse metrics, i.e., ∂ξᵢ/∂xᵢ (or ξₓ). To do this
  # the jacobian matrix has to be assembled and then inverted
  _inverse_metrics_3d(scheme, domain)

  # Next we find the **conservative** metrics which are a pain in the
  # arse and need to be found a la ξ̂x = (y_η z)_ζ − (y_ζ z)_η; 
  # we need mixed derivatives. In this example, we already have 
  # y_η and z, so to find the derivative
  # with respect to ζ, we need to interpolate to k±½, so we can then
  # compute ∂ϕ/∂ζ = (ϕₖ₊½ - ϕₖ₋½) / Δζ.

  # "*" is the operator here, since we're taking the ∂(yη * z)/∂ζ and
  # we already have yη and z in memory. If it was ∂(yη + z)/∂ζ, the 
  # operator would be "+", and so on.

  # ξ̂x = (y_η z)_ζ − (y_ζ z)_η
  _iterp_mixed_terms_to_edge(scheme, y_η, zc, y_ηz_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _iterp_mixed_terms_to_edge(scheme, y_ζ, zc, y_ζz_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η 
  _conserved_metric_term!(ξ̂x, (y_ηz_ζₖ₊½, y_ζz_ηⱼ₊½), (ζ_ax, η_ax), domain) # -> ξ̂x

  # η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
  _iterp_mixed_terms_to_edge(scheme, y_ζ, zc, y_ζz_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _iterp_mixed_terms_to_edge(scheme, y_ξ, zc, y_ξz_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _conserved_metric_term!(η̂x, (y_ζz_ξᵢ₊½, y_ξz_ζₖ₊½), (ξ_ax, ζ_ax), domain) # -> η̂x

  # ζ̂x = (y_ξ z)_η − (y_η z)_ξ
  _iterp_mixed_terms_to_edge(scheme, y_ξ, zc, y_ξz_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η
  _iterp_mixed_terms_to_edge(scheme, y_η, zc, y_ηz_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _conserved_metric_term!(ζ̂x, (y_ξz_ηⱼ₊½, y_ηz_ξᵢ₊½), (η_ax, ξ_ax), domain) # -> ζ̂x

  # ξ̂y = (z_η x)_ζ − (z_ζ x)_η
  _iterp_mixed_terms_to_edge(scheme, z_η, xc, z_ηx_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _iterp_mixed_terms_to_edge(scheme, z_ζ, xc, z_ζx_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η 
  _conserved_metric_term!(ξ̂y, (z_ηx_ζₖ₊½, z_ζx_ηⱼ₊½), (ζ_ax, η_ax), domain) # -> ξ̂y

  # η̂y = (z_ζ x)_ξ − (z_ξ x)_ζ
  _iterp_mixed_terms_to_edge(scheme, z_ζ, xc, z_ζx_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _iterp_mixed_terms_to_edge(scheme, z_ξ, xc, z_ξx_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _conserved_metric_term!(η̂y, (z_ζx_ξᵢ₊½, z_ξx_ζₖ₊½), (ξ_ax, ζ_ax), domain) # -> η̂y

  # ζ̂y = (z_ξ x)_η − (z_η x)_ξ
  _iterp_mixed_terms_to_edge(scheme, z_ξ, xc, z_ξx_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η
  _iterp_mixed_terms_to_edge(scheme, z_η, xc, z_ηx_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _conserved_metric_term!(ζ̂y, (z_ξx_ηⱼ₊½, z_ηx_ξᵢ₊½), (η_ax, ξ_ax), domain) # -> ζ̂y

  # ξ̂z = (x_η y)_ζ − (x_ζ y)_η
  _iterp_mixed_terms_to_edge(scheme, x_η, yc, x_ηy_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _iterp_mixed_terms_to_edge(scheme, x_ζ, yc, x_ζy_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η 
  _conserved_metric_term!(ξ̂z, (x_ηy_ζₖ₊½, x_ζy_ηⱼ₊½), (ζ_ax, η_ax), domain) # -> ξ̂z

  # η̂z = (x_ζ y)_ξ − (x_ξ y)_ζ 
  _iterp_mixed_terms_to_edge(scheme, x_ζ, yc, x_ζy_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _iterp_mixed_terms_to_edge(scheme, x_ξ, yc, x_ξy_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _conserved_metric_term!(η̂z, (x_ζy_ξᵢ₊½, x_ξy_ζₖ₊½), (ξ_ax, ζ_ax), domain) # -> η̂z

  # ζ̂z = (x_ξ y)_η − (x_η y)_ξ
  _iterp_mixed_terms_to_edge(scheme, x_ξ, yc, x_ξy_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η
  _iterp_mixed_terms_to_edge(scheme, x_η, yc, x_ηy_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _conserved_metric_term!(ζ̂z, (x_ξy_ηⱼ₊½, x_ηy_ξᵢ₊½), (η_ax, ξ_ax), domain) # -> ζ̂z

  # Now interpolate the conserved metrics to the edges
  _iterp_to_edge!(scheme, ξ̂x, ξ̂x_i₊½, ξ_ax, domain) # -> ξ̂x_i₊½
  _iterp_to_edge!(scheme, η̂x, η̂x_i₊½, ξ_ax, domain) # -> η̂x_i₊½
  _iterp_to_edge!(scheme, ζ̂x, ζ̂x_i₊½, ξ_ax, domain) # -> ζ̂x_i₊½
  _iterp_to_edge!(scheme, ξ̂y, ξ̂y_i₊½, ξ_ax, domain) # -> ξ̂y_i₊½
  _iterp_to_edge!(scheme, η̂y, η̂y_i₊½, ξ_ax, domain) # -> η̂y_i₊½
  _iterp_to_edge!(scheme, ζ̂y, ζ̂y_i₊½, ξ_ax, domain) # -> ζ̂y_i₊½
  _iterp_to_edge!(scheme, ξ̂z, ξ̂z_i₊½, ξ_ax, domain) # -> ξ̂z_i₊½
  _iterp_to_edge!(scheme, η̂z, η̂z_i₊½, ξ_ax, domain) # -> η̂z_i₊½
  _iterp_to_edge!(scheme, ζ̂z, ζ̂z_i₊½, ξ_ax, domain) # -> ζ̂z_i₊½

  _iterp_to_edge!(scheme, ξ̂x, ξ̂x_j₊½, η_ax, domain) # -> ξ̂x_j₊½
  _iterp_to_edge!(scheme, η̂x, η̂x_j₊½, η_ax, domain) # -> η̂x_j₊½
  _iterp_to_edge!(scheme, ζ̂x, ζ̂x_j₊½, η_ax, domain) # -> ζ̂x_j₊½
  _iterp_to_edge!(scheme, ξ̂y, ξ̂y_j₊½, η_ax, domain) # -> ξ̂y_j₊½
  _iterp_to_edge!(scheme, η̂y, η̂y_j₊½, η_ax, domain) # -> η̂y_j₊½
  _iterp_to_edge!(scheme, ζ̂y, ζ̂y_j₊½, η_ax, domain) # -> ζ̂y_j₊½
  _iterp_to_edge!(scheme, ξ̂z, ξ̂z_j₊½, η_ax, domain) # -> ξ̂z_j₊½
  _iterp_to_edge!(scheme, η̂z, η̂z_j₊½, η_ax, domain) # -> η̂z_j₊½
  _iterp_to_edge!(scheme, ζ̂z, ζ̂z_j₊½, η_ax, domain) # -> ζ̂z_j₊½

  _iterp_to_edge!(scheme, ξ̂x, ξ̂x_k₊½, ζ_ax, domain) # -> ξ̂x_k₊½ 
  _iterp_to_edge!(scheme, η̂x, η̂x_k₊½, ζ_ax, domain) # -> η̂x_k₊½ 
  _iterp_to_edge!(scheme, ζ̂x, ζ̂x_k₊½, ζ_ax, domain) # -> ζ̂x_k₊½ 
  _iterp_to_edge!(scheme, ξ̂y, ξ̂y_k₊½, ζ_ax, domain) # -> ξ̂y_k₊½ 
  _iterp_to_edge!(scheme, η̂y, η̂y_k₊½, ζ_ax, domain) # -> η̂y_k₊½ 
  _iterp_to_edge!(scheme, ζ̂y, ζ̂y_k₊½, ζ_ax, domain) # -> ζ̂y_k₊½ 
  _iterp_to_edge!(scheme, ξ̂z, ξ̂z_k₊½, ζ_ax, domain) # -> ξ̂z_k₊½ 
  _iterp_to_edge!(scheme, η̂z, η̂z_k₊½, ζ_ax, domain) # -> η̂z_k₊½ 
  _iterp_to_edge!(scheme, ζ̂z, ζ̂z_k₊½, ζ_ax, domain) # -> ζ̂z_k₊½ 

  return nothing
end

# 2D version
"""
    update_metrics!(scheme, xc, yc, domain) 

Update the 2D grid metrics given the cell-centroid positions `xc` and `yc`
"""
function update_metrics!(
  scheme::MEG6Scheme{2}, xc::AbstractArray{T,2}, yc::AbstractArray{T,2}, domain
) where {T}
  # compute the ∂x/∂ξ|ᵢ, ∂x/∂η|ᵢ terms

  # Each dimension is on a different axes in the various metric arrays.
  # These are used since the metric functions are defined for
  # arbitrary dimensions, so we need to define what axis we're working on
  ξ_ax, η_ax = (1, 2)

  _cell_center_metric!(
    scheme.∂x∂ξᵢ.∂x₁.∂ξ₁, xc, scheme.∂_cache, scheme.∂²_cache, ξ_ax, domain
  ) # ∂x/∂ξ
  _cell_center_metric!(
    scheme.∂x∂ξᵢ.∂x₁.∂ξ₂, xc, scheme.∂_cache, scheme.∂²_cache, η_ax, domain
  ) # ∂x/∂η
  _cell_center_metric!(
    scheme.∂x∂ξᵢ.∂x₂.∂ξ₁, yc, scheme.∂_cache, scheme.∂²_cache, ξ_ax, domain
  ) # ∂y/∂ξ
  _cell_center_metric!(
    scheme.∂x∂ξᵢ.∂x₂.∂ξ₂, yc, scheme.∂_cache, scheme.∂²_cache, η_ax, domain
  ) # ∂y/∂η

  # compute the inverse metrics, i.e., ∂ξᵢ/∂xᵢ (or ξₓ). To do this
  # the jacobian matrix has to be assembled and then inverted
  _inverse_metrics(scheme, domain) # -> this also populates scheme.J

  # Next we find the **conservative** metrics. This is much more
  # complicated in 3D...

  # ξ̂x = ξx / J
  ξx = scheme.∂ξ∂xᵢ.∂ξ₁.∂x₁
  ξ̂xᵢ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₁.i₊½
  ξ̂xⱼ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₁.∂x₁.j₊½
  _iterp_mixed_terms_to_edge(
    ξ̂xᵢ₊½, ξx, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, ξ_ax, domain, /
  )
  _iterp_mixed_terms_to_edge(
    ξ̂xⱼ₊½, ξx, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, η_ax, domain, /
  )

  # η̂x = ηx / J
  ηx = scheme.∂ξ∂xᵢ.∂ξ₂.∂x₁
  η̂xᵢ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₁.i₊½
  η̂xⱼ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₁.j₊½
  _iterp_mixed_terms_to_edge(
    η̂xᵢ₊½, ηx, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, ξ_ax, domain, /
  )
  _iterp_mixed_terms_to_edge(
    η̂xⱼ₊½, ηx, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, η_ax, domain, /
  )

  # ξ̂y = ξy / J
  ξy = scheme.∂ξ∂xᵢ.∂ξ₁.∂x₂
  ξ̂yᵢ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₁.i₊½
  ξ̂yⱼ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₁.j₊½
  _iterp_mixed_terms_to_edge(
    ξ̂yᵢ₊½, ξy, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, ξ_ax, domain, /
  )
  _iterp_mixed_terms_to_edge(
    ξ̂yⱼ₊½, ξy, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, η_ax, domain, /
  )

  # η̂y = ηy / J
  ηy = scheme.∂ξ∂xᵢ.∂ξ₂.∂x₂
  η̂yᵢ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₂.i₊½
  η̂yⱼ₊½ = scheme.∂ξ̂∂xᵢ₊½.∂ξ₂.∂x₂.j₊½
  _iterp_mixed_terms_to_edge(
    η̂yᵢ₊½, ηy, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, ξ_ax, domain, /
  )
  _iterp_mixed_terms_to_edge(
    η̂yⱼ₊½, ηy, scheme.Jᵢ, scheme.∂_cache, scheme.∂²_cache, η_ax, domain, /
  )

  return nothing
end

# 1D version
function update_metrics!(scheme::MEG6Scheme{1}, xc::AbstractArray{T,1}, domain) where {T}
  # compute the ∂x/∂ξ|ᵢ, ∂x/∂η|ᵢ terms

  # Each dimension is on a different axes in the various metric arrays.
  # These are used since the metric functions are defined for
  # arbitrary dimensions, so we need to define what axis we're working on
  ξ_ax = 1

  _cell_center_metric(scheme, xc, scheme.∂xᵢ∂ξᵢ.x.ξ, ξ_ax, domain) # ∂x/∂ξ

  # compute the inverse metrics, i.e., ∂ξᵢ/∂xᵢ (or ξₓ). To do this
  # the jacobian matrix has to be assembled and then inverted
  _inverse_metrics_1d(scheme, domain)

  # Next we find the **conservative** metrics 

  # ξ̂x = ξx / J
  _iterp_mixed_terms_to_edge(scheme, ξx, J, ξ̂xᵢ₊½, ξ_ax, domain, /)

  return nothing
end

# ----------------------------------------------------------------
# Private functions not intended to be used outside of this module
# ----------------------------------------------------------------

# 1st derivative operator
@inline function forward_∂_2nd_order(ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂ϕ = -(3 / 2) * ϕ[0] + 2ϕ[+1] - (1 / 2) * ϕ[+2]
  return _∂ϕ
end

@inline function backward_∂_2nd_order(ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂ϕ = +(3 / 2) * ϕ[0] - 2ϕ[-1] + (1 / 2) * ϕ[-2]
  return _∂ϕ
end

@inline function forward_mixed_∂_4th_order(ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂ϕ = (1 / 12) * (-25ϕ[0] + 48ϕ[+1] - 36ϕ[+2] + 16ϕ[+3] - 3ϕ[+4])
  # _∂ϕ = (
  #   -(1 / 4) * ϕ[0] - (5 / 6) * ϕ[+1] + (3 / 2) * ϕ[+2] - (1 / 2) * ϕ[+3] + (1 / 12) * ϕ[+4]
  # )
  return _∂ϕ
end

@inline function backward_mixed_∂_4th_order(ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂ϕ = (1 / 12) * (+25ϕ[0] - 48ϕ[-1] + 36ϕ[-2] - 16ϕ[-3] + 3ϕ[-4])
  # _∂ϕ = (
  #   -(1 / 4) * ϕ[0] - (5 / 6) * ϕ[-1] + (3 / 2) * ϕ[-2] - (1 / 2) * ϕ[-3] + (1 / 12) * ϕ[-4]
  # )
  return _∂ϕ
end

@inline function central_∂_4th_order(ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂ϕ = ((2 / 3) * (ϕ[+1] - ϕ[-1]) - (1 / 12) * (ϕ[+2] - ϕ[-2]))
  return _∂ϕ
end

@inline function ∂(ϕ::OffsetVector{T,SVector{7,T}}) where {T}
  f1 = 3 / 4
  f2 = 3 / 20
  f3 = 1 / 60
  _∂ϕ = f1 * (ϕ[+1] - ϕ[-1]) - f2 * (ϕ[+2] - ϕ[-2]) + f3 * (ϕ[+3] - ϕ[-3])
  return _∂ϕ
end

# 2nd derivative operator
@inline function ∂²(
  ϕ::OffsetVector{T,SVector{3,T}}, ∂ϕ::OffsetVector{T,SVector{3,T}}
) where {T}
  _∂²ϕ = 2(ϕ[+1] - 2ϕ[0] + ϕ[-1]) - (∂ϕ[+1] - ∂ϕ[-1]) / 2
  return _∂²ϕ
end

@inline function forward_∂²(∂ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂²ϕ = -(3 / 2) * ∂ϕ[+0] + 2∂ϕ[+1] - (1 / 2) * ∂ϕ[+2]
  return _∂²ϕ
end

@inline function forward_mixed_∂²(∂ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂²ϕ =
    -(1 / 4) * ∂ϕ[0] - (5 / 6) * ∂ϕ[+1] + (3 / 2) * ∂ϕ[+2] - (1 / 2) * ∂ϕ[+3] +
    (1 / 12) * ∂ϕ[+4]
  return _∂²ϕ
end

@inline function backward_∂²(∂ϕ::OffsetVector{T,SVector{3,T}}) where {T}
  _∂²ϕ = +(3 / 2) * ∂ϕ[+0] - 2∂ϕ[-1] + (1 / 2) * ∂ϕ[-2]
  return _∂²ϕ
end

@inline function backward_mixed_∂²(∂ϕ::OffsetVector{T,SVector{5,T}}) where {T}
  _∂²ϕ =
    +(1 / 4) * ∂ϕ[0] + (5 / 6) * ∂ϕ[-1] - (3 / 2) * ∂ϕ[-2] + (1 / 2) * ∂ϕ[-3] -
    (1 / 12) * ∂ϕ[-4]
  return _∂²ϕ
end

# @inline function ∂²_central(
#   ϕ::OffsetVector{T,SVector{3,T}}, ∂ϕ::OffsetVector{T,SVector{3,T}}
# ) where {T}
#   # _∂²ϕ = -(1 / 2) * (∂ϕ[+1] - ∂ϕ[-1]) - 2(ϕ[+1] - 2ϕ[0] + ϕ[-1])
#   _∂²ϕ = -(1 / 2) * (∂ϕ[+1] - ∂ϕ[-1]) - 2(ϕ[+1] - 2ϕ[0] + ϕ[-1])
#   return _∂²ϕ
# end

function _get_inner_gradient(
  ϕ::AbstractArray{T,N}, ∂ϕ∂ξ::AbstractArray{T,N}, axis::Int, domain
) where {T,N}
  return nothing
end

function _boundary_∂!(
  ∂ϕ::AbstractArray{T,N}, ϕ::AbstractArray{T,N}, axis::Int, domain
) where {T,N}

  # Create CartesianIndices iterators for the thin boundary regions 
  # along the given axis. For example, lets say the boundary axis is "i",
  # or axis=1, so the lower_b1 corresponts to CartesianIndices((1:1, 1:10))
  # for a domain of CartesianIndices((1:10,1:10)), lower_b2 is (2:2, 1:10) 
  # and so on. The same logic applies to upper_b1, etc.
  lower_b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis
  lower_b2 = lower_boundary_indices(domain, axis, +1) # first index + 1 on given boundary axis
  lower_b3 = lower_boundary_indices(domain, axis, +2) # first index + 2 on given boundary axis
  upper_b3 = upper_boundary_indices(domain, axis, -2) # last index on given boundary axis
  upper_b2 = upper_boundary_indices(domain, axis, -1) # last index + 1 on given boundary axis
  upper_b1 = upper_boundary_indices(domain, axis, 0)  # last index + 2 on given boundary axis

  # Do 1st derivatives first...
  # one-sided forward ∂ϕ
  for idx in lower_b1
    imp = (0, 2) # i.e offset of 0:2 
    offset_idx = plus_minus(idx, axis, imp)
    _ϕ = OffsetVector(SVector{3}(view(ϕ, offset_idx)), 0:2)
    ∂ϕ[idx] = forward_∂_2nd_order(_ϕ)
  end

  for idx in lower_b2
    imp = (0, 4) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    _ϕ = OffsetVector(SVector{5}(view(ϕ, offset_idx)), 0:4)
    ∂ϕ[idx] = forward_mixed_∂_4th_order(_ϕ)
  end

  for idx in lower_b3
    imp = (2, 2) # i.e offset of -2:2
    offset_idx = plus_minus(idx, axis, imp)
    _ϕ = OffsetVector(SVector{5}(view(ϕ, offset_idx)), -2:2)
    ∂ϕ[idx] = central_∂_4th_order(_ϕ)
  end

  # one-sided backward ∂ϕ
  for idx in upper_b3
    imp = (2, 2) # i.e offset of -2:2
    offset_idx = plus_minus(idx, axis, imp)
    _ϕ = OffsetVector(SVector{5}(view(ϕ, offset_idx)), -2:2)
    ∂ϕ[idx] = central_∂_4th_order(_ϕ)
  end

  for idx in upper_b2
    imp = (4, 0) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    _ϕ = OffsetVector(SVector{5}(view(ϕ, offset_idx)), -4:0)
    ∂ϕ[idx] = backward_mixed_∂_4th_order(_ϕ)
  end

  for idx in upper_b1
    imp = (2, 0) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    _ϕ = OffsetVector(SVector{3}(view(ϕ, offset_idx)), -2:0)
    ∂ϕ[idx] = backward_∂_2nd_order(_ϕ)
  end

  return nothing
end

function _boundary_∂AB!(
  ∂AB::AbstractArray{T,N},
  A::AbstractArray{T,N},
  B::AbstractArray{T,N},
  axis::Int,
  domain,
  operator=*,
) where {T,N}

  # Create CartesianIndices iterators for the thin boundary regions 
  # along the given axis. For example, lets say the boundary axis is "i",
  # or axis=1, so the lower_b1 corresponts to CartesianIndices((1:1, 1:10))
  # for a domain of CartesianIndices((1:10,1:10)), lower_b2 is (2:2, 1:10) 
  # and so on. The same logic applies to upper_b1, etc.
  lower_b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis
  lower_b2 = lower_boundary_indices(domain, axis, +1) # first index + 1 on given boundary axis
  lower_b3 = lower_boundary_indices(domain, axis, +2) # first index + 2 on given boundary axis
  upper_b3 = upper_boundary_indices(domain, axis, -2) # last index on given boundary axis
  upper_b2 = upper_boundary_indices(domain, axis, -1) # last index + 1 on given boundary axis
  upper_b1 = upper_boundary_indices(domain, axis, 0)  # last index + 2 on given boundary axis

  # Do 1st derivatives first...
  # one-sided forward ∂AB
  for idx in lower_b1
    imp = (0, 2) # i.e offset of 0:2 
    offset_idx = plus_minus(idx, axis, imp)
    _A = SVector{3}(view(A, offset_idx))
    _B = SVector{3}(view(B, offset_idx))
    _AB = OffsetVector(operator.(_A, _B), 0:2)
    ∂AB[idx] = forward_∂_2nd_order(_AB)
  end

  for idx in lower_b2
    imp = (0, 4) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    _A = SVector{5}(view(A, offset_idx))
    _B = SVector{5}(view(B, offset_idx))
    _AB = OffsetVector(operator.(_A, _B), 0:4)
    ∂AB[idx] = forward_mixed_∂_4th_order(_AB)
  end

  for idx in lower_b3
    imp = (2, 2) # i.e offset of -2:2
    offset_idx = plus_minus(idx, axis, imp)
    _A = SVector{5}(view(A, offset_idx))
    _B = SVector{5}(view(B, offset_idx))
    _AB = OffsetVector(operator.(_A, _B), -2:2)
    ∂AB[idx] = central_∂_4th_order(_AB)
  end

  # one-sided backward ∂AB
  for idx in upper_b3
    imp = (2, 2) # i.e offset of -2:2
    offset_idx = plus_minus(idx, axis, imp)
    # _ϕ = OffsetVector(SVector{5}(view(ϕ, offset_idx)), -2:2)
    _A = SVector{5}(view(A, offset_idx))
    _B = SVector{5}(view(B, offset_idx))
    _AB = OffsetVector(operator.(_A, _B), -2:2)
    ∂AB[idx] = central_∂_4th_order(_AB)
  end

  for idx in upper_b2
    imp = (4, 0) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    # _ϕ = OffsetVector(SVector{5}(view(ϕ, offset_idx)), -4:0)
    _A = SVector{5}(view(A, offset_idx))
    _B = SVector{5}(view(B, offset_idx))
    _AB = OffsetVector(operator.(_A, _B), -4:0)
    ∂AB[idx] = backward_mixed_∂_4th_order(_AB)
  end

  for idx in upper_b1
    imp = (2, 0) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    # _ϕ = OffsetVector(SVector{3}(view(ϕ, offset_idx)), -2:0)
    _A = SVector{3}(view(A, offset_idx))
    _B = SVector{3}(view(B, offset_idx))
    _AB = OffsetVector(operator.(_A, _B), -2:0)
    ∂AB[idx] = backward_∂_2nd_order(_AB)
  end

  return nothing
end

function _boundary_∂²!(
  ∂²ϕ::AbstractArray{T,N}, ∂ϕ::AbstractArray{T,N}, axis::Int, domain
) where {T,N}

  # Create CartesianIndices iterators for the thin boundary regions 
  # along the given axis. For example, lets say the boundary axis is "i",
  # or axis=1, so the lower_b1 corresponts to CartesianIndices((1:1, 1:10))
  # for a domain of CartesianIndices((1:10,1:10)), lower_b2 is (2:2, 1:10) 
  # and so on. The same logic applies to upper_b1, etc.
  lower_b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis
  lower_b2 = lower_boundary_indices(domain, axis, +1) # first index + 1 on given boundary axis
  # lower_b3 = lower_boundary_indices(domain, axis, +2) # first index + 2 on given boundary axis
  # upper_b3 = upper_boundary_indices(domain, axis, -2) # last index on given boundary axis
  upper_b2 = upper_boundary_indices(domain, axis, -1) # last index + 1 on given boundary axis
  upper_b1 = upper_boundary_indices(domain, axis, 0)  # last index + 2 on given boundary axis

  # Now do 2nd derivatives, since they depend on the 1st derivatives
  # one-sided forward ∂²ϕ
  for idx in lower_b1
    imp = (0, 2) # i.e offset of 0:2
    offset_idx = plus_minus(idx, axis, imp)
    _∂ϕ = OffsetVector(SVector{3}(view(∂ϕ, offset_idx)), 0:2)
    ∂²ϕ[idx] = forward_∂²(_∂ϕ)
  end

  for idx in lower_b2
    imp = (0, 4) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    _∂ϕ = OffsetVector(SVector{5}(view(∂ϕ, offset_idx)), 0:4)
    ∂²ϕ[idx] = forward_mixed_∂²(_∂ϕ)
  end

  # for idx in lower_b3
  #   imp = (1, 1) # i.e offset of -2:2
  #   offset_idx = plus_minus(idx, axis, imp)
  #   _ϕ = OffsetVector(SVector{3}(view(ϕ, offset_idx)), -1:1)
  #   _∂ϕ = OffsetVector(SVector{3}(view(∂ϕ, offset_idx)), -1:1)
  #   ∂²ϕ[idx] = ∂²_central(_ϕ, _∂ϕ)
  # end

  # one-sided backward ∂²ϕ
  # for idx in upper_b3
  #   imp = (1, 1) # i.e offset of -1:1
  #   offset_idx = plus_minus(idx, axis, imp)
  #   _ϕ = OffsetVector(SVector{3}(view(ϕ, offset_idx)), -1:1)
  #   _∂ϕ = OffsetVector(SVector{3}(view(∂ϕ, offset_idx)), -1:1)
  #   ∂²ϕ[idx] = ∂²_central(_ϕ, _∂ϕ)
  # end

  for idx in upper_b2
    imp = (4, 0) # i.e offset of 0:4
    offset_idx = plus_minus(idx, axis, imp)
    _∂ϕ = OffsetVector(SVector{5}(view(∂ϕ, offset_idx)), -4:0)
    ∂²ϕ[idx] = backward_mixed_∂²(_∂ϕ)
  end

  for idx in upper_b1
    imp = (2, 0) # i.e offset of 0:2
    offset_idx = plus_minus(idx, axis, imp)
    _∂ϕ = OffsetVector(SVector{3}(view(∂ϕ, offset_idx)), -2:0)
    ∂²ϕ[idx] = backward_∂²(_∂ϕ)
  end

  return nothing
end

@inline function _conserved_metric_term!(
  ∂ϕ, (ϕ1, ϕ2)::NTuple{2,AbstractArray{T,3}}, (ax1, ax2)::NTuple{2,Int}, domain
) where {T}

  # All the conservative metrics look like the following
  # for a 3D mesh. 

  # ξ̂x = (y_η z)_ζ − (y_ζ z)_η
  # η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
  # ζ̂x = (y_ξ z)_η − (y_η z)_ξ

  # ξ̂y = (z_η x)_ζ − (z_ζ x)_η
  # η̂y = (z_ζ x)_ξ − (z_ξ x)_ζ
  # ζ̂y = (z_ξ x)_η − (z_η x)_ξ

  # ξ̂z = (x_η y)_ζ − (x_ζ y)_η
  # η̂z = (x_ζ y)_ξ − (x_ξ y)_ζ
  # ζ̂z = (x_ξ y)_η − (x_η y)_ξ

  # This function takes the edge based arrays
  # and does the differencing to find the outer derivative.
  # Taking ξ̂x = (y_η z)_ζ − (y_ζ z)_η as an example, y_η z
  # is defined at k+1/2 for all cells, and to get (y_η z)_ζ
  # we do (y_η z)_ζ = ((y_η z)ₖ₊½ - (y_η z)ₖ₋½) / Δζ

  # In the loop below, ϕ1 is (y_η z)ₖ₊½ and ϕ2 is (y_ζ z)ⱼ₊½

  # This is a stupidly simple loop, but it works for arbitrary dimensions
  @inbounds for idx in domain1
    # Sticking with the example for (y_η z)_ζ = ((y_η z)ₖ₊½ - (y_η z)ₖ₋½) / Δζ,
    # but this will work for arbitrary axes
    k₋½ = down(idx, ax1, 1)# i.e. k₋½ is simply k₊½ at [i, j, k-1]
    j₋½ = down(idx, ax2, 1)
    j₊½ = k₊½ = idx

    ∂ϕ[idx] = (
      (ϕ1[k₊½] - ϕ1[k₋½]) - # / Δζ is in the formula, but Δξ = Δη = Δζ = 1 by definition
      (ϕ2[j₊½] - ϕ2[j₋½])   # / Δη is in the formula, but Δξ = Δη = Δζ = 1 by definition
    )
    # TODO: do we need an epsilon check here like ∂ϕ[idx] = ∂ϕ[idx] * abs(∂ϕ[idx] >= eps) ?
  end

  return nothing
end

# 3D version
function _inverse_metrics(scheme::MEG6Scheme{3}, domain)
  @inbounds for idx in domain
    xξ = scheme.∂x∂ξᵢ.x₁.ξ₁[idx] # ∂x/∂ξ
    yξ = scheme.∂x∂ξᵢ.x₂.ξ₁[idx] # ∂y/∂ξ
    zξ = scheme.∂x∂ξᵢ.x₃.ξ₁[idx] # ∂z/∂ξ
    xη = scheme.∂x∂ξᵢ.x₁.ξ₂[idx] # ∂x/∂η
    yη = scheme.∂x∂ξᵢ.x₂.ξ₂[idx] # ∂y/∂η
    zη = scheme.∂x∂ξᵢ.x₃.ξ₂[idx] # ∂z/∂η
    xζ = scheme.∂x∂ξᵢ.x₁.ξ₃[idx] # ∂x/∂ζ
    yζ = scheme.∂x∂ξᵢ.x₂.ξ₃[idx] # ∂y/∂ζ
    zζ = scheme.∂x∂ξᵢ.x₃.ξ₃[idx] # ∂z/∂ζ

    # inverse jacobian matrix
    J⁻¹ = @SMatrix [
      xξ yξ zξ
      xη yη zη
      xζ yζ zζ
    ]

    J = inv(J⁻¹) # jacobian matrix
    scheme.Jᵢ[idx] = det(J) # "the Jacobian"... why can't we use different names??

    scheme.∂ξ∂xᵢ.ξ₁.x₁[idx] = J[1, 1] # ∂ξ/∂x
    scheme.∂ξ∂xᵢ.ξ₁.x₂[idx] = J[1, 2] # ∂ξ/∂y
    scheme.∂ξ∂xᵢ.ξ₁.x₃[idx] = J[1, 3] # ∂ξ/∂z
    scheme.∂ξ∂xᵢ.ξ₂.x₁[idx] = J[2, 1] # ∂η/∂x
    scheme.∂ξ∂xᵢ.ξ₂.x₂[idx] = J[2, 2] # ∂η/∂y
    scheme.∂ξ∂xᵢ.ξ₂.x₃[idx] = J[2, 3] # ∂η/∂z
    scheme.∂ξ∂xᵢ.ξ₃.x₁[idx] = J[3, 1] # ∂ζ/∂x
    scheme.∂ξ∂xᵢ.ξ₃.x₂[idx] = J[3, 2] # ∂ζ/∂y
    scheme.∂ξ∂xᵢ.ξ₃.x₃[idx] = J[3, 3] # ∂ζ/∂z
  end

  return nothing
end

# 2D version
function _inverse_metrics(scheme::MEG6Scheme{2}, domain)
  # extended_domain = expand(domain, 1)
  @inbounds for idx in domain
    xξ = scheme.∂x∂ξᵢ.∂x₁.∂ξ₁[idx] # ∂x/∂ξ
    xη = scheme.∂x∂ξᵢ.∂x₁.∂ξ₂[idx] # ∂x/∂η
    yξ = scheme.∂x∂ξᵢ.∂x₂.∂ξ₁[idx] # ∂y/∂ξ
    yη = scheme.∂x∂ξᵢ.∂x₂.∂ξ₂[idx] # ∂y/∂η

    # inverse jacobian matrix
    J⁻¹ = @SMatrix [
      xξ yξ
      xη yη
    ]

    J = inv(J⁻¹)
    scheme.Jᵢ[idx] = det(J) # "the Jacobian".... why can't we use different names??

    scheme.∂ξ∂xᵢ.∂ξ₁.∂x₁[idx] = J[1, 1] # ∂ξ/∂x
    scheme.∂ξ∂xᵢ.∂ξ₁.∂x₂[idx] = J[1, 2] # ∂ξ/∂y
    scheme.∂ξ∂xᵢ.∂ξ₂.∂x₁[idx] = J[2, 1] # ∂η/∂x
    scheme.∂ξ∂xᵢ.∂ξ₂.∂x₂[idx] = J[2, 2] # ∂η/∂y
  end

  return nothing
end

# 1D version
function _inverse_metrics(scheme::MEG6Scheme{1}, domain)
  # 1D is extremely simple! No need for matrices, inv(), or det()
  @inbounds for idx in domain
    J⁻¹ = scheme.∂x∂ξᵢ.x₁.ξ₁[idx] # ∂x/∂ξ
    J = 1 / J⁻¹
    scheme.Jᵢ[idx] = J # "the Jacobian".... why can't we use different names??
    scheme.∂ξ∂xᵢ.ξ₁.x₁[idx] = J # ∂ξ/∂x
  end

  return nothing
end

"""
Interpolate cell-centered combination of A and B terms to the edge. Consistency 
requires that the whole term is interpolated, i.e, A*B, rather than interpolating 
them separately and then applying the operator.
"""
function _iterp_mixed_terms_to_edge(
  ABᵢ₊½::AbstractArray{T,N},
  A::AbstractArray{T,N},
  B::AbstractArray{T,N},
  ∂AB, # cache arrays
  ∂²AB, # cache arrays
  axis::Int, # 1 = i, 2 = j, 3 = k
  domain, # domain iterator, e.g., CartesianIndices
  operator=*, # we're interpolating operator.(A,B), e.g, A*B
) where {T,N}

  # intermediate-arrays
  # ∂AB = scheme.∂_cache
  # ∂²AB = scheme.∂²_cache
  # TODO: fix this shit!
  # gradients for the boundaries require mixed and one-sided stencils
  # _boundary_∂AB!(∂AB, A, B, axis, domain)

  # since we already did the special boundary treatment, 
  # shrink the domain where we apply the standard stencils below
  inner_domain_m3 = expand(domain, axis, -3) # shrink by 3 boundary cells

  # find the derivative term ∂AB, or ∂(A*B)/∂i
  for idx in inner_domain_m3
    # get i-3:i+3 on whichever axis, e.g. i, j, or k
    offset = 3
    offset_idx = plus_minus(idx, axis, offset)

    # create small stencil vectors
    A_stencil = SVector{2offset + 1}(view(A, offset_idx))
    B_stencil = SVector{2offset + 1}(view(B, offset_idx))

    # The mixed term is usually A*B, with the * operator (which we apply element-wise with the .)
    # This could be any 2-argument function, i.e. +,-,*,/, etc.
    # Use an offset vector so the cell in question is always at 0,
    # and its convienient to use -n:n indices for stencil operations
    AB_stencil = OffsetVector(operator.(A_stencil, B_stencil), (-offset):offset)
    ∂AB[idx] = ∂(AB_stencil)
  end

  # gradients for the boundaries require mixed and one-sided stencils
  # _boundary_∂²!(∂²AB, ∂AB, axis, domain)

  # cell-centered ∂²ϕ
  inner_domain_m1 = expand(domain, axis, -1) # shrink by 3 boundary cells

  # find the 2nd derivative term ∂²AB, or ∂²(A*B)/∂i
  for idx in inner_domain_m1
    offset = 1
    offset_idx = plus_minus(idx, axis, offset)

    A_stencil = SVector{3}(view(A, offset_idx))
    B_stencil = SVector{3}(view(B, offset_idx))
    AB_stencil = OffsetVector(operator.(A_stencil, B_stencil), (-offset):offset)
    ∂AB_stencil = OffsetVector(SVector{3}(view(∂AB, offset_idx)), (-offset):offset)
    ∂²AB[idx] = ∂²(AB_stencil, ∂AB_stencil)
  end

  # # treat boundaries separately for the edge terms
  # lower_b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis
  # for idx in lower_b1
  #   up_one = up(idx, axis, 1) # i.e. [i+1, j, k]
  #   down_one = down(idx, axis, 1) # i.e. [i+1, j, k]

  #   AB = operator.(A[idx], B[idx]) # i.e. A[i,j,k] * B[i,j,k]
  #   AB_up_one = operator.(A[up_one], B[up_one]) # i.e, A[i+1, j, k] * B[i+1, j, k]

  #   ABᴸᵢ₊½ = AB + 0.5∂AB[idx] + ∂²AB[idx] / 12
  #   ABᴿᵢ₊½ = AB_up_one - 0.5∂AB[up_one] + ∂²AB[up_one] / 12

  #   ABᵢ₋½ = AB - 0.5∂AB[idx] + ∂²AB[idx] / 12

  #   ABᵢ₊½[down_one] = ABᵢ₋½ # we need to make sure that ᵢ₊½ is stored along the low boundary
  #   ABᵢ₊½[idx] = 0.5(ABᴸᵢ₊½ + ABᴿᵢ₊½)
  # end

  # upper_b1 = upper_boundary_indices(domain, axis, 0)  # last index on given boundary axis
  # for idx in upper_b1
  #   # down_one = down(idx, axis, 1) # i.e. [i-1, j, k]

  #   AB = operator.(A[idx], B[idx]) # i.e. A[i,j,k] * B[i,j,k]
  #   # AB_down_one = operator.(A[down_one], B[down_one]) # i.e, A[i+1, j, k] * B[i+1, j, k]

  #   ABᴸᵢ₊½ = AB + 0.5∂AB[idx] + ∂²AB[idx] / 12

  #   # ABᴸᵢ₋½ = AB_down_one + 0.5∂AB[down_one] + ∂²AB[down_one] / 12
  #   # ABᴿᵢ₋½ = AB - 0.5∂AB[idx] + ∂²AB[idx] / 12

  #   ABᵢ₊½[idx] = ABᴸᵢ₊½
  # end

  # interpolate to the edge
  for idx in domain
    up_one = up(idx, axis, 1)# i.e. [i+1, j, k]
    AB = operator.(A[idx], B[idx]) # i.e. A[i,j,k] * B[i,j,k]
    AB_up_one = operator.(A[up_one], B[up_one]) # i.e, A[i+1, j, k] * B[i+1, j, k]

    # the name says ᵢ₊½, but it could be i₊½, j₊½, or k₊½, it just means the "+" edge
    ABᴸᵢ₊½ = AB + 0.5∂AB[idx] + ∂²AB[idx] / 12
    ABᴿᵢ₊½ = AB_up_one - 0.5∂AB[up_one] + ∂²AB[up_one] / 12
    ABᵢ₊½[idx] = 0.5(ABᴸᵢ₊½ + ABᴿᵢ₊½)
  end

  return nothing
end

"""
Interpolate cell-centered quantity `A` to cell edges for arbitary dimensions
"""
function _iterp_to_edge!(
  scheme::MEG6Scheme{N,T},
  A::AbstractArray{T,N},
  Aᵢ₊½::AbstractArray{T,N},
  axis::Int, # 1 = i, 2 = j, 3 = k
  domain, # domain iterator, e.g., CartesianIndices
) where {T,N}

  # intermediate-arrays
  ∂A = scheme.∂_cache
  ∂²A = scheme.∂²_cache

  # find the derivative term ∂A
  @inbounds for idx in domain
    # get i-3:i+3 on whichever axis, e.g. i, j, or k
    offset = 3
    offset_idx = plus_minus(idx, axis, offset)

    # create small stencil vectors
    A_stencil = OffsetVector(
      SVector{2offset + 1}(view(A, offset_idx)), # stencil
      (-offset):offset,
    )

    # Use an offset vector so the cell in question is always at 0,
    # and its convienient to use -n:n indices for stencil operatios
    ∂A[idx] = ∂ϕ(A_stencil)
  end

  # find the 2nd derivative term ∂²A
  @inbounds for idx in domain
    offset = 1
    offset_idx = plus_minus(idx, axis, offset)

    A_stencil = OffsetVector(
      SVector{2offset + 1}(view(A, offset_idx)), # stencil  
      (-offset):offset,
    )

    ∂A_stencil = OffsetVector(
      SVector{2offset + 1}(view(∂A, offset_idx)), # stencil  
      (-offset):offset,
    )

    ∂²A[idx] = ∂²ϕ(A_stencil, ∂A_stencil)
  end

  # interpolate to the edge using a central scheme (no upwind biasing)
  # so we simply average the L and R terms
  @inbounds for idx in domain
    up_one = up(idx, axis, 1)# i.e. [i+1, j, k]

    # the name says ᵢ₊½, but it could be i₊½, j₊½, or k₊½, it just means the "+" edge
    Aᴸᵢ₊½ = A[idx] + 0.5∂A[idx] + ∂²A[idx] / 12
    Aᴿᵢ₊½ = A[up_one] - 0.5∂A[up_one] + ∂²A[up_one] / 12
    Aᵢ₊½[idx] = 0.5(Aᴸᵢ₊½ + Aᴿᵢ₊½)
  end

  return nothing
end

"""
Compute a cell-centered term ∂ϕ∂ξ based on the MEG6 gradient scheme. Provide the axis to determine
which direction
"""
function _cell_center_metric!(
  ∂ϕ∂ξ, ϕ::AbstractArray{T,N}, ∂ϕ, ∂²ϕ, axis::Int, domain
) where {T,N}
  # intermediate-arrays
  # ∂ϕ = scheme.∂_cache
  # ∂²ϕ = scheme.∂²_cache

  # gradients for the boundaries require mixed and one-sided stencils
  _boundary_∂!(∂ϕ, ϕ, axis, domain)

  # since we already did the special boundary treatment, 
  # shrink the domain where we apply the standard stencils below
  inner_domain_m3 = expand(domain, axis, -3) # shrink by 3 boundary cells

  # cell-centered ∂ϕ
  for idx in inner_domain_m3
    offset = 3
    # get i-3:i+3 on whichever axis, e.g. i, j, or k
    offset_idx = plus_minus(idx, axis, offset)

    # Use an small offset vector so the cell in question is always at 0,
    # and its convienient to use -n:n indices for stencil operatios
    _ϕ = OffsetVector(
      SVector{2offset + 1}(view(ϕ, offset_idx)), # stencil for [i-3:i+3,j,k] or whichever axis
      (-offset):offset,
    )
    ∂ϕ[idx] = ∂(_ϕ)
  end

  # gradients for the boundaries require mixed and one-sided stencils
  _boundary_∂²!(∂²ϕ, ∂ϕ, axis, domain)

  # cell-centered ∂²ϕ
  inner_domain_m1 = expand(domain, axis, -1) # shrink by 3 boundary cells

  for idx in inner_domain_m1
    offset = 1
    # get i-1:i+1 on whichever axis, e.g. i, j, or k
    offset_idx = plus_minus(idx, axis, offset)

    _ϕ = OffsetVector(
      SVector{2offset + 1}(view(ϕ, offset_idx)), # stencil
      (-offset):offset,
    )
    _∂ϕ = OffsetVector(
      SVector{2offset + 1}(view(∂ϕ, offset_idx)), # stencil
      (-offset):offset,
    )
    ∂²ϕ[idx] = ∂²(_ϕ, _∂ϕ)
  end

  # ------------------------------------------------------------------
  # now find ∂ϕ∂ξ
  # ------------------------------------------------------------------
  # the edge names say ᵢ₊½, but it could be i₊½, j₊½, or k₊½, it just means the "+/-" edge

  # treat boundaries separately for the edge terms
  lower_b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis
  for idx in lower_b1
    up_one = up(idx, axis, 1) # i.e. [i+1, j, k]

    ϕᴸᵢ₊½ = ϕ[idx] + 0.5∂ϕ[idx] + ∂²ϕ[idx] / 12
    ϕᴿᵢ₊½ = ϕ[up_one] - 0.5∂ϕ[up_one] + ∂²ϕ[up_one] / 12

    ϕᵢ₋½ = ϕ[idx] - 0.5∂ϕ[idx] + ∂²ϕ[idx] / 12

    ϕᵢ₊½ = 0.5(ϕᴸᵢ₊½ + ϕᴿᵢ₊½)

    ∂ϕ∂ξ[idx] = ϕᵢ₊½ - ϕᵢ₋½ # / Δξ, but Δξ = 1 by definition
  end

  upper_b1 = upper_boundary_indices(domain, axis, 0)  # last index on given boundary axis
  for idx in upper_b1
    down_one = down(idx, axis, 1) # i.e. [i-1, j, k]

    ϕᵢ₊½ = ϕ[idx] + 0.5∂ϕ[idx] + ∂²ϕ[idx] / 12

    ϕᴸᵢ₋½ = ϕ[down_one] + 0.5∂ϕ[down_one] + ∂²ϕ[down_one] / 12
    ϕᴿᵢ₋½ = ϕ[idx] - 0.5∂ϕ[idx] + ∂²ϕ[idx] / 12

    ϕᵢ₋½ = 0.5(ϕᴸᵢ₋½ + ϕᴿᵢ₋½)

    ∂ϕ∂ξ[idx] = ϕᵢ₊½ - ϕᵢ₋½ # / Δξ, but Δξ = 1 by definition
  end

  inner_domain_m1 = expand(domain, axis, -1) # shrink by 1 boundary cell
  for idx in inner_domain_m1
    up_one = up(idx, axis, 1) # i.e. [i+1, j, k]
    down_one = down(idx, axis, 1) # i.e. [i-1, j, k]

    ϕᴸᵢ₊½ = ϕ[idx] + 0.5∂ϕ[idx] + ∂²ϕ[idx] / 12
    ϕᴸᵢ₋½ = ϕ[down_one] + 0.5∂ϕ[down_one] + ∂²ϕ[down_one] / 12

    ϕᴿᵢ₊½ = ϕ[up_one] - 0.5∂ϕ[up_one] + ∂²ϕ[up_one] / 12
    ϕᴿᵢ₋½ = ϕ[idx] - 0.5∂ϕ[idx] + ∂²ϕ[idx] / 12

    ϕᵢ₊½ = 0.5(ϕᴸᵢ₊½ + ϕᴿᵢ₊½)
    ϕᵢ₋½ = 0.5(ϕᴸᵢ₋½ + ϕᴿᵢ₋½)

    ∂ϕ∂ξ[idx] = ϕᵢ₊½ - ϕᵢ₋½ # / Δξ, but Δξ = 1 by definition
  end

  return nothing
end
