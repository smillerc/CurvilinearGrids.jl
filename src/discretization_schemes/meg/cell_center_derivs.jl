"""
    cell_center_derivatives!(scheme::MonotoneExplicitGradientScheme, ϕ_ξ, ϕ, axis, inner_cell_domain; compute_gradients=false, use_one_sided_on_edges=true, ϵ=50eps(eltype(ϕ)))
    cell_center_derivatives!(scheme::MonotoneExplicitGradientScheme, ϕ_ξ, ∂²ϕ, ∂ϕ, ϕ, axis, inner_cell_domain; compute_gradients=false, use_one_sided_on_edges=true, ϵ=50eps(eltype(ϕ)))

Compute the monotone explicit gradient approximation to the cell-centered derivative `ϕ_ξ` of the quantity `ϕ` along the provided `axis`. When cache arrays `∂ϕ` and `∂²ϕ` are supplied, they are reused; otherwise cached buffers on the scheme are used internally.

Set `compute_gradients=true` to populate the derivative buffers before reconstructing the edge values, and use `use_one_sided_on_edges=false` to force central differences everywhere. The `ϵ` keyword zeroes derivatives with magnitude below the tolerance.
"""
function cell_center_derivatives!(
  scheme::MonotoneExplicitGradientScheme,
  ϕ_ξ::AbstractArray{T,N},
  ϕ::AbstractArray{T,N},
  axis::Int,
  inner_cell_domain;
  kwargs...,
) where {T,N}
  cell_center_derivatives!(
    scheme, ϕ_ξ, scheme.cache.∂²ϕ, scheme.cache.∂ϕ, ϕ, axis, inner_cell_domain; kwargs...
  )
end

function cell_center_derivatives!(
  scheme::MonotoneExplicitGradientScheme,
  ϕ_ξ::AbstractArray{T,N},
  ∂²ϕ::AbstractArray{T,N},
  ∂ϕ::AbstractArray{T,N},
  ϕ::AbstractArray{T,N},
  axis::Int,
  inner_cell_domain;
  compute_gradients=false,
  use_one_sided_on_edges=true,
  ϵ=50eps(eltype(ϕ)),
) where {T,N}
  backend = KernelAbstractions.get_backend(scheme.cache.outer_deriv_1)

  if compute_gradients
    compute_first_derivatives!(
      scheme, ∂ϕ, ϕ, axis, inner_cell_domain, backend, use_one_sided_on_edges
    )
    compute_second_derivatives!(
      scheme, ∂²ϕ, ∂ϕ, ϕ, axis, inner_cell_domain, backend, use_one_sided_on_edges
    )
  end

  if use_one_sided_on_edges
    # shrink the derivative domain so the central deriv form can be used
    inner_domain = expand(inner_cell_domain, axis, -1)
  else
    # otherwise, use the central deriv form everywhere
    inner_domain = inner_cell_domain
  end

  meg_inner_kernel!(backend)(
    ϕ, ∂ϕ, ∂²ϕ, ϕ_ξ, ϵ, inner_domain, axis; ndrange=size(inner_domain)
  )

  if use_one_sided_on_edges
    index_offset = 0

    # select first index along given boundary axis
    lo_domain = lower_boundary_indices(inner_cell_domain, axis, index_offset)

    # the lo-side derivative ∂ϕ/∂ξ can only use ϕᴿᵢ₋½ instead of ϕᵢ₋½
    meg_lo_edge_kernel!(backend)(
      ϕ, ∂ϕ, ∂²ϕ, ϕ_ξ, ϵ, lo_domain, axis; ndrange=size(lo_domain)
    )

    # select last index along given boundary axis
    hi_domain = upper_boundary_indices(inner_cell_domain, axis, index_offset)

    # the hi-side derivative ∂ϕ/∂ξ can only use ϕᴸᵢ₊½ instead of ϕᵢ₊½
    meg_hi_edge_kernel!(backend)(
      ϕ, ∂ϕ, ∂²ϕ, ϕ_ξ, ϵ, hi_domain, axis; ndrange=size(hi_domain)
    )
  end

  return nothing
end

@kernel inbounds = true function meg_inner_kernel!(ϕ, ∂ϕ, ∂²ϕ, ϕ_ξ, ϵ, domain, axis)
  idx = @index(Global, Linear)
  i = domain[idx]

  ᵢ₋₁ = shift(i, axis, -1)
  ᵢ₊₁ = shift(i, axis, +1)

  # average the L/R reconstructed values
  ϕᵢ₊½ = (ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i]) + ϕᴿᵢ₊½(ϕ[ᵢ₊₁], ∂ϕ[ᵢ₊₁], ∂²ϕ[ᵢ₊₁])) / 2
  ϕᵢ₋½ = (ϕᴸᵢ₊½(ϕ[ᵢ₋₁], ∂ϕ[ᵢ₋₁], ∂²ϕ[ᵢ₋₁]) + ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])) / 2

  _ϕ_ξ = ϕᵢ₊½ - ϕᵢ₋½
  _ϕ_ξ = _ϕ_ξ * (abs(_ϕ_ξ) >= ϵ)
  ϕ_ξ[i] = _ϕ_ξ
end

@kernel inbounds = true function meg_lo_edge_kernel!(ϕ, ∂ϕ, ∂²ϕ, ϕ_ξ, ϵ, domain, axis)
  idx = @index(Global, Linear)
  i = domain[idx]

  ᵢ₊₁ = shift(i, axis, +1)

  ϕᵢ₊½ = (ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i]) + ϕᴿᵢ₊½(ϕ[ᵢ₊₁], ∂ϕ[ᵢ₊₁], ∂²ϕ[ᵢ₊₁])) / 2
  ϕᵢ₋½ = ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])

  _ϕ_ξ = ϕᵢ₊½ - ϕᵢ₋½
  _ϕ_ξ = _ϕ_ξ * (abs(_ϕ_ξ) >= ϵ)
  ϕ_ξ[i] = _ϕ_ξ
end

@kernel inbounds = true function meg_hi_edge_kernel!(ϕ, ∂ϕ, ∂²ϕ, ϕ_ξ, ϵ, domain, axis)
  idx = @index(Global, Linear)
  i = domain[idx]

  ᵢ₋₁ = shift(i, axis, -1)

  ϕᵢ₊½ = ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])
  ϕᵢ₋½ = (ϕᴸᵢ₊½(ϕ[ᵢ₋₁], ∂ϕ[ᵢ₋₁], ∂²ϕ[ᵢ₋₁]) + ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])) / 2

  _ϕ_ξ = ϕᵢ₊½ - ϕᵢ₋½
  _ϕ_ξ = _ϕ_ξ * (abs(_ϕ_ξ) >= ϵ)
  ϕ_ξ[i] = _ϕ_ξ
end
