"""
    interpolate_to_edges!(scheme::MonotoneExplicitGradientScheme, ϕᵢ₊½::AbstractArray{T,N}, ϕ::AbstractArray{T,N}, axis::Int, domain, nhalo=3) where {T,N}

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
  scheme::MonotoneExplicitGradientScheme,
  ϕᵢ₊½::AbstractArray{T,N},
  ϕ::AbstractArray{T,N},
  axis::Int,
  domain,
  backend,
) where {T,N}

  #
  if CartesianIndices(ϕᵢ₊½) != CartesianIndices(ϕ)
    error("The given ϕᵢ₊½ and ϕ arrays must have the same size and indexing scheme")
  end

  ∂ϕ = scheme.cache.∂ϕ
  ∂²ϕ = scheme.cache.∂²ϕ
  compute_first_derivatives!(scheme, ∂ϕ, ϕ, axis, domain, backend)
  compute_second_derivatives!(scheme, ∂²ϕ, scheme.cache.∂ϕ, ϕ, axis, domain, backend)

  if size(domain) == size(ϕ)
    inner_domain = expand(domain, axis, -1)
    lo_domain = lower_boundary_indices(domain, axis, +1) # select first index along given boundary axis
    hi_domain = upper_boundary_indices(domain, axis, -1)  # select last index along given boundary axis
  else
    inner_domain = domain
    lo_domain = lower_boundary_indices(domain, axis, 0)
    hi_domain = upper_boundary_indices(domain, axis, 0)
  end

  inner_kernel!(backend)(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, inner_domain, axis; ndrange=size(inner_domain))
  lo_edge_kernel!(backend)(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, lo_domain, axis; ndrange=size(lo_domain))
  hi_edge_kernel!(backend)(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, hi_domain, axis; ndrange=size(hi_domain))

  return nothing
end

@kernel inbounds = true function inner_kernel!(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, domain, axis)
  idx = @index(Global, Linear)
  i = domain[idx]

  ᵢ₊₁ = shift(i, axis, +1)

  # average the L/R reconstructed values
  ϕᵢ₊½[i] = (
    ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i]) +     # left state
    ϕᴿᵢ₊½(ϕ[ᵢ₊₁], ∂ϕ[ᵢ₊₁], ∂²ϕ[ᵢ₊₁]) # right state
  ) / 2
end

@kernel inbounds = true function lo_edge_kernel!(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, lo_domain, axis)
  idx = @index(Global, Linear)
  i = lo_domain[idx]

  # the lo-side derivative ∂ϕ/∂ξ can only use ϕᴿᵢ₋½ instead of ϕᵢ₋½
  ᵢ₋₁ = shift(i, axis, -1)
  ϕᵢ₊½[ᵢ₋₁] = ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])
end

@kernel inbounds = true function hi_edge_kernel!(ϕᵢ₊½, ϕ, ∂ϕ, ∂²ϕ, hi_domain, axis)
  idx = @index(Global, Linear)
  i = hi_domain[idx]

  # the hi-side derivative ∂ϕ/∂ξ can only use ϕᴸᵢ₊½ instead of ϕᵢ₊½
  ϕᵢ₊½[i] = ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])
end
