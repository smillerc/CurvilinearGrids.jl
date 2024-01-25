using OffsetArrays
using StaticArrays

dim_names(::NTuple{1,Int}) = (:ξ,)
dim_names(::NTuple{2,Int}) = (:ξ, :η)
dim_names(::NTuple{3,Int}) = (:ξ, :η, :ζ)

coord_names(::NTuple{1,Int}) = (:x₁,)
coord_names(::NTuple{2,Int}) = (:x₁, :x₂)
coord_names(::NTuple{3,Int}) = (:x₁, :x₂, :x₃)

function metric_vars(dims)
  vars = coord_names(dims)
  gradients = NamedTuple{(vars)}(
    cell_centered_vars(dims) for i in 1:length(vars) # (:x=(ξ,η,ζ), :y=(ξ,η,ζ), ...)
  )

  return gradients
end

function inv_metric_vars(dims)
  vars = coord_names(dims)
  gradients = NamedTuple{(vars)}(
    cell_centered_vars(dims) for i in 1:length(vars) # (:x=(ξ,η,ζ), :y=(ξ,η,ζ), ...)
  )

  return gradients
end

function cell_centered_vars(dims)
  ξηζ = dim_names(dims)
  # (:ξ=Array, :η=Array, :ζ=Array)
  return NamedTuple{(ξηζ)}(zeros(dims) for i in 1:length(ξηζ))
end

struct MEG6Scheme{A,B,C,D}
  ∂ϕᵢ::A # gradient terms of any scalar term ϕ (these are cached and re-used for each metric term)
  ∂²ϕᵢ::B # gradient terms of any scalar term ϕ (these are cached and re-used for each metric term)
  ∂ξᵢ∂xᵢ::C # metrics; ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  ∂ξ̂ᵢ∂xᵢ::C # conservative metrics, where ξ̂ = ξ/J; ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  ∂ξ̂ᵢ∂xᵢ_edge::C # conservative metrics at the cell edges ((i,j,k)+1/2), where ξ̂ = ξ/J; ∂(ξ,η,ζ)/∂z, ∂(ξ,η,ζ)∂y, ∂(ξ,η,ζ)∂z
  ∂xᵢ∂ξᵢ::D # ∂x/∂(ξ,η,ζ), ∂y/∂(ξ,η,ζ), ∂z/∂(ξ,η,ζ)
  J
end

function MEG6Scheme(celldims::NTuple{N,Int}) where {N}
  ∂ϕᵢ = cell_centered_vars(celldims)
  ∂²ϕᵢ = cell_centered_vars(celldims)
  x_ξ = metric_vars(dims)
  ξ_x = inv_metric_vars(dims)
  return MEG6Scheme(∂ϕᵢ, ∂²ϕᵢ, x_ξ, ξ_x)
end

dims(m::MEG6Scheme) = keys(m.∂ϕᵢ)

@inline function ∂ϕ(ϕ::OffsetVector{T,SVector{7,T}}) where {T}
  f1 = 3 / 4
  f2 = 3 / 20
  f3 = 1 / 60
  _∂ϕ = f1 * (ϕ[+1] - ϕ[-1]) - f2 * (ϕ[+2] - ϕ[-2]) + f3 * (ϕ[+3] - ϕ[-3])
  return _∂ϕ
end

# ∂²ϕ(ϕ, ∂ϕ) = 2(ϕ[+1] - 2ϕ[0] + ϕ[-1]) - f1 * (∂ϕ[+1] - ∂ϕ[-1])
@inline function ∂²ϕ(
  ϕ::OffsetVector{T,SVector{3,T}}, ∂ϕ::OffsetVector{T,SVector{3,T}}
) where {T}
  _∂²ϕ = 2(ϕ[+1] - 2ϕ[0] + ϕ[-1]) - (∂ϕ[+1] - ∂ϕ[-1]) / 2
  return _∂²ϕ
end

function update_metrics!(scheme::MEG6Scheme, x, y, z)
  # compute the ∂x/∂ξ|ᵢ, ∂x/∂η|ᵢ terms

  # Each dimension is on a different axes in the various metric arrays.
  # These are used since the metric functions are defined for
  # arbitrary dimensions, so we need to define what axis we're working on
  ξ_ax, η_ax, ζ_ax = (1, 2, 3)

  _cell_center_metric(scheme, x, scheme.∂xᵢ∂ξᵢ.x.ξ, ξ_ax) # ∂x/∂ξ
  _cell_center_metric(scheme, x, scheme.∂xᵢ∂ξᵢ.x.η, η_ax) # ∂x/∂η
  _cell_center_metric(scheme, x, scheme.∂xᵢ∂ξᵢ.x.ζ, ζ_ax) # ∂x/∂ζ

  _cell_center_metric(scheme, y, scheme.∂xᵢ∂ξᵢ.y.ξ, ξ_ax) # ∂y/∂ξ
  _cell_center_metric(scheme, y, scheme.∂xᵢ∂ξᵢ.y.η, η_ax) # ∂y/∂η
  _cell_center_metric(scheme, y, scheme.∂xᵢ∂ξᵢ.y.ζ, ζ_ax) # ∂y/∂ζ

  _cell_center_metric(scheme, z, scheme.∂xᵢ∂ξᵢ.z.ξ, ξ_ax) # ∂z/∂ξ
  _cell_center_metric(scheme, z, scheme.∂xᵢ∂ξᵢ.z.η, η_ax) # ∂z/∂η
  _cell_center_metric(scheme, z, scheme.∂xᵢ∂ξᵢ.z.ζ, ζ_ax) # ∂z/∂ζ

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
  _iterp_mixed_terms_to_edge(scheme, y_η, z, y_ηz_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _iterp_mixed_terms_to_edge(scheme, y_ζ, z, y_ζz_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η 
  _conserved_metric_term!(ξ̂x, (y_ηz_ζₖ₊½, y_ζz_ηⱼ₊½), (ζ_ax, η_ax), domain) # -> ξ̂x

  # η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
  _iterp_mixed_terms_to_edge(scheme, y_ζ, z, y_ζz_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _iterp_mixed_terms_to_edge(scheme, y_ξ, z, y_ξz_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _conserved_metric_term!(η̂x, (y_ζz_ξᵢ₊½, y_ξz_ζₖ₊½), (ξ_ax, ζ_ax), domain) # -> η̂x

  # ζ̂x = (y_ξ z)_η − (y_η z)_ξ
  _iterp_mixed_terms_to_edge(scheme, y_ξ, z, y_ξz_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η
  _iterp_mixed_terms_to_edge(scheme, y_η, z, y_ηz_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _conserved_metric_term!(ζ̂x, (y_ξz_ηⱼ₊½, y_ηz_ξᵢ₊½), (η_ax, ξ_ax), domain) # -> ζ̂x

  # ξ̂y = (z_η x)_ζ − (z_ζ x)_η
  _iterp_mixed_terms_to_edge(scheme, z_η, x, z_ηx_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _iterp_mixed_terms_to_edge(scheme, z_ζ, x, z_ζx_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η 
  _conserved_metric_term!(ξ̂y, (z_ηx_ζₖ₊½, z_ζx_ηⱼ₊½), (ζ_ax, η_ax), domain) # -> ξ̂y

  # η̂y = (z_ζ x)_ξ − (z_ξ x)_ζ
  _iterp_mixed_terms_to_edge(scheme, z_ζ, x, z_ζx_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _iterp_mixed_terms_to_edge(scheme, z_ξ, x, z_ξx_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _conserved_metric_term!(η̂y, (z_ζx_ξᵢ₊½, z_ξx_ζₖ₊½), (ξ_ax, ζ_ax), domain) # -> η̂y

  # ζ̂y = (z_ξ x)_η − (z_η x)_ξ
  _iterp_mixed_terms_to_edge(scheme, z_ξ, x, z_ξx_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η
  _iterp_mixed_terms_to_edge(scheme, z_η, x, z_ηx_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _conserved_metric_term!(ζ̂y, (z_ξx_ηⱼ₊½, z_ηx_ξᵢ₊½), (η_ax, ξ_ax), domain) # -> ζ̂y

  # ξ̂z = (x_η y)_ζ − (x_ζ y)_η
  _iterp_mixed_terms_to_edge(scheme, x_η, y, x_ηy_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _iterp_mixed_terms_to_edge(scheme, x_ζ, y, x_ζy_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η 
  _conserved_metric_term!(ξ̂z, (x_ηy_ζₖ₊½, x_ζy_ηⱼ₊½), (ζ_ax, η_ax), domain) # -> ξ̂z

  # η̂z = (x_ζ y)_ξ − (x_ξ y)_ζ 
  _iterp_mixed_terms_to_edge(scheme, x_ζ, y, x_ζy_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _iterp_mixed_terms_to_edge(scheme, x_ξ, y, x_ξy_ζₖ₊½, ζ_ax, domain, *) # interp to k+½, since the outer deriv is in ζ
  _conserved_metric_term!(η̂z, (x_ζy_ξᵢ₊½, x_ξy_ζₖ₊½), (ξ_ax, ζ_ax), domain) # -> η̂z

  # ζ̂z = (x_ξ y)_η − (x_η y)_ξ
  _iterp_mixed_terms_to_edge(scheme, x_ξ, y, x_ξy_ηⱼ₊½, η_ax, domain, *) # interp to j+½, since the outer deriv is in η
  _iterp_mixed_terms_to_edge(scheme, x_η, y, x_ηy_ξᵢ₊½, ξ_ax, domain, *) # interp to i+½, since the outer deriv is in ξ
  _conserved_metric_term!(ζ̂z, (x_ξy_ηⱼ₊½, x_ηy_ξᵢ₊½), (η_ax, ξ_ax), domain) # -> ζ̂z

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

# 2D version
function _inverse_metrics_2d(scheme, domain)
  @inbounds for idx in domain
    xξ = scheme.∂xᵢ∂ξᵢ.x.ξ[idx] # ∂x/∂ξ
    yξ = scheme.∂xᵢ∂ξᵢ.y.ξ[idx] # ∂y/∂ξ
    xη = scheme.∂xᵢ∂ξᵢ.x.η[idx] # ∂x/∂η
    yη = scheme.∂xᵢ∂ξᵢ.y.η[idx] # ∂y/∂η

    # inverse jacobian matrix
    J⁻¹ = @SMatrix [
      xξ yξ
      xη yη
    ]

    J = inv(J⁻¹)
    scheme.J[idx] = det(J) # "the Jacobian".... why can't we use different names??

    scheme.∂ξᵢ∂xᵢ.ξ.x₁[idx] = J[1, 1] # ∂ξ/∂x
    scheme.∂ξᵢ∂xᵢ.ξ.x₂[idx] = J[1, 2] # ∂ξ/∂y
    scheme.∂ξᵢ∂xᵢ.η.x₁[idx] = J[2, 1] # ∂η/∂x
    scheme.∂ξᵢ∂xᵢ.η.x₂[idx] = J[2, 2] # ∂η/∂y
  end

  return nothing
end

# 3D version
function _inverse_metrics_3d(scheme, domain)
  @inbounds for idx in domain
    xξ = scheme.∂xᵢ∂ξᵢ.x.ξ[idx] # ∂x/∂ξ
    yξ = scheme.∂xᵢ∂ξᵢ.y.ξ[idx] # ∂y/∂ξ
    zξ = scheme.∂xᵢ∂ξᵢ.z.ξ[idx] # ∂z/∂ξ
    xη = scheme.∂xᵢ∂ξᵢ.x.η[idx] # ∂x/∂η
    yη = scheme.∂xᵢ∂ξᵢ.y.η[idx] # ∂y/∂η
    zη = scheme.∂xᵢ∂ξᵢ.z.η[idx] # ∂z/∂η
    xζ = scheme.∂xᵢ∂ξᵢ.x.ζ[idx] # ∂x/∂ζ
    yζ = scheme.∂xᵢ∂ξᵢ.y.ζ[idx] # ∂y/∂ζ
    zζ = scheme.∂xᵢ∂ξᵢ.z.ζ[idx] # ∂z/∂ζ

    # inverse jacobian matrix
    J⁻¹ = @SMatrix [
      xξ yξ zξ
      xη yη zη
      xζ yζ zζ
    ]

    J = inv(J⁻¹) # jacobian matrix
    scheme.J[idx] = det(J) # "the Jacobian"... why can't we use different names??

    scheme.∂ξᵢ∂xᵢ.ξ.x₁[idx] = J[1, 1] # ∂ξ/∂x
    scheme.∂ξᵢ∂xᵢ.ξ.x₂[idx] = J[1, 2] # ∂ξ/∂y
    scheme.∂ξᵢ∂xᵢ.ξ.x₃[idx] = J[1, 3] # ∂ξ/∂z
    scheme.∂ξᵢ∂xᵢ.η.x₁[idx] = J[2, 1] # ∂η/∂x
    scheme.∂ξᵢ∂xᵢ.η.x₂[idx] = J[2, 2] # ∂η/∂y
    scheme.∂ξᵢ∂xᵢ.η.x₃[idx] = J[2, 3] # ∂η/∂z
    scheme.∂ξᵢ∂xᵢ.ζ.x₁[idx] = J[3, 1] # ∂ζ/∂x
    scheme.∂ξᵢ∂xᵢ.ζ.x₂[idx] = J[3, 2] # ∂ζ/∂y
    scheme.∂ξᵢ∂xᵢ.ζ.x₃[idx] = J[3, 3] # ∂ζ/∂z
  end

  return nothing
end

# function _interp_conserv_metrics_to_edges(
#   scheme, ϕ::AbstractArray{T,2}, J::AbstractArray{T,2}, ∂ξ̂∂xᵢ₊½, ∂ξ̂∂xⱼ₊½
# ) where {T}

#   # interpolate from the cell-centered quantity ϕ to the conservative
#   # edge metric term ϕ/J|ᵢ₊½, where  ϕ is whatever metric quantity, like ξₓ
#   # ϕ̂ is ϕ / J, e.g. the conservative metric, like ξ̂ₓ = ξₓ / J

#   ∂ϕ̂ᵢ = scheme.∂ϕᵢ.ξ
#   ∂ϕ̂ⱼ = scheme.∂ϕᵢ.η
#   ∂²ϕ̂ᵢ = scheme.∂²ϕᵢ.ξ
#   ∂²ϕ̂ⱼ = scheme.∂²ϕᵢ.η

#   # edge metrics we're updating
#   ∂ϕ̂∂xᵢ₊½ = ∂ξᵢ∂ϕ̂ᵢ.x
#   ∂ϕ̂∂yⱼ₊½ = ∂ξᵢ∂ϕ̂ᵢ.y

#   # find the derivative term ∂ϕ̂ᵢ, or ∂(ϕ/J)/∂i
#   @inbounds for idx in mesh.iterators.cell.domain
#     i, j = idx.I
#     ϕ_stencil = SVector{7}(view(ϕ, (i - 3):(i + 3), j))
#     J_stencil = SVector{7}(view(J, (i - 3):(i + 3), j))
#     ϕ̂ = OffsetVector(ϕ_stencil / J_stencil, -3:3)
#     ∂ϕ̂ᵢ[idx] = ∂ϕ(ϕ̂)
#   end

#   @inbounds for idx in mesh.iterators.cell.domain
#     i, j = idx.I
#     ϕ_stencil = SVector{7}(view(ϕ, i, (j - 3):(j + 3)))
#     J_stencil = SVector{7}(view(J, i, (j - 3):(j + 3)))
#     ϕ̂ = OffsetVector(ϕ_stencil / J_stencil, -3:3)
#     ∂ϕ̂ⱼ[idx] = ∂ϕ(ϕ̂)
#   end

#   @inbounds for idx in mesh.iterators.cell.domain
#     i, j = idx.I

#     ϕ_stencil = SVector{3}(view(ϕ, (i - 1):(i + 1), j))
#     J_stencil = SVector{3}(view(J, (i - 1):(i + 1), j))
#     ϕ̂ = OffsetVector(ϕ_stencil / J_stencil, -1:1)
#     ∂ϕ̂ = OffsetVector(SVector{3}(view(∂ϕ̂ⱼ, (i - 1):(i + 1), j)), -1:1)
#     ∂²ϕ̂ᵢ[idx] = ∂²ϕ(ϕ̂, ∂ϕ̂)
#   end

#   @inbounds for idx in mesh.iterators.cell.domain
#     i, j = idx.I
#     ϕ_stencil = SVector{3}(view(ϕ, i, (j - 1):(j + 1)))
#     J_stencil = SVector{3}(view(J, i, (j - 1):(j + 1)))
#     ϕ̂ = OffsetVector(ϕ_stencil / J_stencil, -1:1)
#     ∂ϕ̂ = OffsetVector(SVector{3}(view(∂ϕ̂ⱼ, i, (j - 1):(j + 1))), -1:1)
#     ∂²ϕ̂ⱼ[idx] = ∂²ϕ(ϕ̂, ∂ϕ̂)
#   end

#   # interpolate to the edges
#   @inbounds for idx in mesh.iterators.cell.domain
#     i, j = idx.I

#     #TODO: how to handle edge ∂ terms?

#     ϕ̂ᴸᵢ₊½ = (ϕ[i, j] / J[i, j]) + 0.5∂ϕ̂ᵢ[i, j] + ∂²ϕ̂ᵢ[i, j] / 12
#     ϕ̂ᴸⱼ₊½ = (ϕ[i, j] / J[i, j]) + 0.5∂ϕ̂ⱼ[i, j] + ∂²ϕ̂ⱼ[i, j] / 12

#     ϕ̂ᴿᵢ₊½ = (ϕ[i + 1, j] / J[i + 1, j]) - 0.5∂ϕ̂ᵢ[i + 1, j] + ∂²ϕ̂ᵢ[i + 1, j] / 12
#     ϕ̂ᴿⱼ₊½ = (ϕ[i, j + 1] / J[i, j + 1]) - 0.5∂ϕ̂ⱼ[i, j + 1] + ∂²ϕ̂ⱼ[i, j + 1] / 12

#     ∂ϕ̂∂xᵢ₊½[idx] = 0.5(ϕ̂ᴸᵢ₊½ + ϕ̂ᴿᵢ₊½)
#     ∂ϕ̂∂yⱼ₊½[idx] = 0.5(ϕ̂ᴸⱼ₊½ + ϕ̂ᴿⱼ₊½)
#   end

#   return nothing
# end

"""
Interpolate cell-centered combination of A and B terms to the edge. Consistency 
requires that the whole term is interpolated, i.e, A*B, rather than interpolating 
them separately and then applying the operator.
"""
function _iterp_mixed_terms_to_edge(
  scheme,
  A::AbstractArray{T,N},
  B::AbstractArray{T,N},
  ABᵢ₊½::AbstractArray{T,N},
  axis::Int, # 1 = i, 2 = j, 3 = k
  domain, # domain iterator, e.g., CartesianIndices
  operator=*, # we're interpolating operator(A,B), e.g, A*B
) where {T}

  # intermediate-arrays
  ∂AB = scheme.∂_cache
  ∂²AB = scheme.∂²_cache

  # find the derivative term ∂AB, or ∂(A*B)/∂i
  @inbounds for idx in domain
    # get i-3:i+3 on whichever axis, e.g. i, j, or k
    offset = 3
    offset_idx = plus_minus(idx, axis, offset)

    # create small stencil vectors
    A_stencil = SVector{2offset + 1}(view(A, offset_idx))
    B_stencil = SVector{2offset + 1}(view(B, offset_idx))

    # The mixed term is usually A*B, with the * operator
    # This could be any 2-argument function, i.e. +,-,*,/, etc.
    # Use an offset vector so the cell in question is always at 0,
    # and its convienient to use -n:n indices for stencil operatios
    AB_stencil = OffsetVector(operator(A_stencil, B_stencil), (-offset):offset)
    ∂AB[idx] = ∂ϕ(AB_stencil)
  end

  # find the 2nd derivative term ∂²AB, or ∂²(A*B)/∂i
  @inbounds for idx in domain
    offset = 1
    offset_idx = plus_minus(idx, axis, offset)

    A_stencil = SVector{3}(view(A, offset_idx))
    B_stencil = SVector{3}(view(B, offset_idx))
    AB_stencil = OffsetVector(operator(A_stencil, B_stencil), (-offset):offset)
    ∂AB_stencil = SVector{3}(view(∂AB, offset_idx))
    ∂²AB[idx] = ∂²ϕ(AB_stencil, ∂AB_stencil)
  end

  # interpolate to the edge
  @inbounds for idx in domain
    up_one = up(idx, axis, 1)# i.e. [i+1, j, k]
    AB = operator(A[idx], B[idx]) # i.e. A[i,j,k] * B[i,j,k]
    AB_up_one = operator(A[up_one], B[up_one]) # i.e, A[i+1, j, k] * B[i+1, j, k]

    # the name says ᵢ₊½, but it could be i₊½, j₊½, or k₊½, it just means the "+" edge
    ABᴸᵢ₊½ = AB + 0.5∂AB[idx] + ∂²AB[idx] / 12
    ABᴿᵢ₊½ = AB_up_one - 0.5∂AB[up_one] + ∂²AB[up_one] / 12
    ABᵢ₊½[idx] = 0.5(ABᴸᵢ₊½ + ABᴿᵢ₊½)
  end

  return nothing
end

function _iterp_to_edge!(
  scheme,
  A::AbstractArray{T,N},
  Aᵢ₊½::AbstractArray{T,N},
  axis::Int, # 1 = i, 2 = j, 3 = k
  domain, # domain iterator, e.g., CartesianIndices
) where {T}

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
Compute a cell-centered metric term, i.e. ∂x∂ξ. Provide the axis to determine
which direction
"""
function _cell_center_metric(scheme, ϕ::AbstractArray{T,N}, ∂ϕ∂ξ, axis::Int) where {T,N}
  # intermediate-arrays
  ∂ϕ = scheme.∂_cache
  ∂²ϕ = scheme.∂²_cache

  # cell-centered ∂ϕ
  @inbounds for idx in domain
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

  # cell-centered ∂²ϕ
  @inbounds for idx in domain
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
  #TODO: how to handle boundary ∂ terms?
  @inbounds for idx in domain
    up_one = up(idx, axis, 1) # i.e. [i+1, j, k]
    down_one = down(idx, axis, 1) # i.e. [i-1, j, k]

    # the name says ᵢ₊½, but it could be i₊½, j₊½, or k₊½, it just means the "+/-" edge
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

# function _cell_centered_mixed_derivative_outer_ζ(scheme)

#   #(y_η⋅z)_ζ, where the _η means ∂/∂η
#   @inbounds for idx in mesh.iterators.cell.domain
#     i, j, k = idx.I

#     yηzₖ₊½ = yηzₖ₋½ = 0
#     yηz_ζ[i, j, k] = yηzₖ₊½ - yηzₖ₋½ # / Δζ, but Δζ = 1 by definition
#   end

#   return nothing
# end

function cons()
  yξz = (
    (2282 / 2880) * (y[i + 1, j, k] - y[i - 1, j, k]) * z[i, j, k] +
    (-546 / 2880) * (y[i + 2, j, k] - y[i - 2, j, k]) * z[i, j, k] +
    (89 / 2880) * (y[i + 3, j, k] - y[i - 3, j, k]) * z[i, j, k] +
    (-3 / 2880) * (y[i + 4, j, k] - y[i - 4, j, k]) * z[i, j, k] +
    (-1 / 2880) * (y[i + 5, j, k] - y[i - 5, j, k]) * z[i, j, k]
  )

  yηz = (
    (2282 / 2880) * (y[i, j + 1, k] - y[i, j - 1, k]) * z[i, j, k] +
    (-546 / 2880) * (y[i, j + 2, k] - y[i, j - 2, k]) * z[i, j, k] +
    (89 / 2880) * (y[i, j + 3, k] - y[i, j - 3, k]) * z[i, j, k] +
    (-3 / 2880) * (y[i, j + 4, k] - y[i, j - 4, k]) * z[i, j, k] +
    (-1 / 2880) * (y[i, j + 5, k] - y[i, j - 5, k]) * z[i, j, k]
  )

  yζz = (
    (2282 / 2880) * (y[i, j, k + 1] - y[i, j, k - 1]) * z[i, j, k] +
    (-546 / 2880) * (y[i, j, k + 2] - y[i, j, k - 2]) * z[i, j, k] +
    (89 / 2880) * (y[i, j, k + 3] - y[i, j, k - 3]) * z[i, j, k] +
    (-3 / 2880) * (y[i, j, k + 4] - y[i, j, k - 4]) * z[i, j, k] +
    (-1 / 2880) * (y[i, j, k + 5] - y[i, j, k - 5]) * z[i, j, k]
  )

  return nothing
end

# for mixed derivatives, like ξ̂x = (y_η z)_ζ − (y_ζ z)_η, we have 
# the y_η terms already computed, and the normal x,y,z terms as well.
# For the outer derivatives, we need to interpolate the product of 
# the inner terms to the edges and take the derivative there.

# ξ̂x = (y_η z)_ζ − (y_ζ z)_η
# η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
# ζ̂x = (y_ξ z)_η − (y_η z)_ξ

# ξ̂y = (z_η x)_ζ − (z_ζ x)_η
# η̂y = (z_ζ x)_ξ − (z_ξ x)_ζ
# ζ̂y = (z_ξ x)_η − (z_η x)_ξ

# ξ̂z = (x_η y)_ζ − (x_ζ y)_η
# η̂z = (x_ζ y)_ξ − (x_ξ y)_ζ
# ζ̂z = (x_ξ y)_η − (x_η y)_ξ

"""
  Apply a delta function to the cartesian index on a specified axis. For 
  example, `δ(3, CartesianIndex(1,2,3))` will give `CartesianIndex(0,0,1)`.
"""
δ(axis, ::CartesianIndex{N}) where {N} = CartesianIndex(ntuple(j -> j == axis ? 1 : 0, N))

"""
"""
function plus_minus(I::CartesianIndex{N}, axis::Int, n::Int) where {N}
  return (I - n * δ(axis, I)):(I + n * δ(axis, I))
end

function up(I::CartesianIndex{N}, axis::Int, n::Int) where {N}
  return I + n * δ(axis, I)
end

function down(I::CartesianIndex{N}, axis::Int, n::Int) where {N}
  return I - n * δ(axis, I)
end