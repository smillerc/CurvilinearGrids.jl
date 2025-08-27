"""
    interpolate_to_edges!(scheme::MontoneExplicitGradientScheme6thOrder, ϕᵢ₊½::AbstractArray{T,N}, ϕ::AbstractArray{T,N}, axis::Int, domain, nhalo=3) where {T,N}

Interpolate the cell-centered quantity `ϕ` to the i+1/2 edge `ϕᵢ₊½` using the MEG6 scheme. 
The `axis` specifies which dimension, e.g. `1 => ϕᵢ₊½, 2 => ϕⱼ₊½, 3 => ϕₖ₊½`. The domain sizes for this
can be a bit tricky due to the staggered nature of the data. The arrays ϕ and ϕᵢ₊½ are the same size, 
but ϕᵢ₊½ only holds the {ᵢ,ⱼ,ₖ}₊½ edge value, so the lowest boundary will not have a value at {ᵢ,ⱼ,ₖ}₋½. 
This is normally fine since the arrays are padded with halo cells along each dimension (nhalo=5 for the 
MEG6 scheme), and ϕᵢ₋½ at the low boundary is accessed with ϕᵢ₊½[ilo-1,j] where ilo is the lowest non-halo index.

The ϕᵢ₊½ and ϕ arrays passed to this function must be the inner domain region + 1 layer of halo cells. So, if
the full domain is CartesianIndices((10, 10)), and the inner domain (with 2 halo cells for this example) is 
CartesianIndices((3:8, 3:8)), the size passed to this function must be CartesianIndices((2:9, 2:9))
"""
function interpolate_to_edge!(
  scheme::MontoneExplicitGradientScheme6thOrder,
  ϕᵢ₊½::AbstractArray{T,N},
  ϕ::AbstractArray{T,N},
  axis::Int,
  domain,
) where {T,N}

  #
  if CartesianIndices(ϕᵢ₊½) != CartesianIndices(ϕ)
    error("The given ϕᵢ₊½ and ϕ arrays must have the same size and indexing scheme")
  end

  backend = KernelAbstractions.get_backend(scheme.cache.outer_deriv_1)
  #   nhalo = scheme.nhalo_for_derivs # 3
  nhalo = 3
  ∂ϕ = scheme.cache.∂ϕ
  ∂²ϕ = scheme.cache.∂²ϕ
  first_deriv!(∂ϕ, ϕ, axis, domain, backend, nhalo)
  second_deriv!(∂²ϕ, scheme.cache.∂ϕ, ϕ, axis, domain, backend, nhalo)

  # domain = CartesianIndices(ϕ)
  inner_domain = domain # expand(domain, axis, 0)
  inner_kernel!(backend)(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, inner_domain, axis; ndrange=size(inner_domain))

  # select first index along given boundary axis
  lo_domain = lower_boundary_indices(domain, axis, 0)
  lo_edge_kernel!(backend)(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, lo_domain, axis; ndrange=size(lo_domain))

  # select last index along given boundary axis
  hi_domain = upper_boundary_indices(domain, axis, 0)
  hi_edge_kernel!(backend)(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, hi_domain, axis; ndrange=size(hi_domain))

  return nothing
end

@kernel inbounds = true function inner_kernel!(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, domain, axis)
  idx = @index(Global)
  i = domain[idx]

  ᵢ₊₁ = shift(i, axis, +1)

  # average the L/R reconstructed values
  ϕᵢ₊½[i] = (
    ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i]) +     # left state
    ϕᴿᵢ₊½(ϕ[ᵢ₊₁], ∂ϕ[ᵢ₊₁], ∂²ϕ[ᵢ₊₁]) # right state
  ) / 2
end

@kernel inbounds = true function lo_edge_kernel!(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, lo_domain, axis)
  idx = @index(Global)
  i = lo_domain[idx]

  # the lo-side derivative ∂ϕ/∂ξ can only use ϕᴿᵢ₋½ instead of ϕᵢ₋½
  ᵢ₋₁ = shift(i, axis, -1)
  ϕᵢ₊½[ᵢ₋₁] = ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])
end

@kernel inbounds = true function hi_edge_kernel!(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, hi_domain, axis)
  idx = @index(Global)
  i = hi_domain[idx]

  # the hi-side derivative ∂ϕ/∂ξ can only use ϕᴸᵢ₊½ instead of ϕᵢ₊½
  ϕᵢ₊½[i] = ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])
end
