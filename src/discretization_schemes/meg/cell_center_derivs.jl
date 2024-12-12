
function cell_center_derivatives!(
  scheme::MontoneExplicitGradientScheme6thOrder, ϕ_ξ, ϕ, axis, domain; nhalo=3
)
  cell_center_derivatives!(
    scheme,
    ϕ_ξ,
    scheme.cache.∂²ϕ,
    scheme.cache.∂ϕ,
    ϕ,
    axis,
    domain;
    compute_gradients=true,
    nhalo=nhalo,
  )
end

function cell_center_derivatives!(
  ::MontoneExplicitGradientScheme6thOrder,
  ϕ_ξ,
  ∂²ϕ,
  ∂ϕ,
  ϕ,
  axis,
  domain;
  compute_gradients=false,
  nhalo=3,
  ϵ=50eps(eltype(ϕ)),
)
  if compute_gradients
    first_deriv!(∂ϕ, ϕ, axis, domain, nhalo)
    second_deriv!(∂²ϕ, ∂ϕ, ϕ, axis, domain, nhalo)
  end

  # # are we computing the derivatives on the inner portion of the domain only? if so
  # # we don't need to use one-sided stencils along the edges
  # if size(domain) == size(ϕ)
  #   inner_domain_only = false
  # domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, axis, -1)
  # inner_domain = domain
  # else
  #   inner_domain_only = true
  #   inner_domain = domain
  # end

  for i in inner_domain
    ᵢ₋₁ = shift(i, axis, -1)
    ᵢ₊₁ = shift(i, axis, +1)

    # average the L/R reconstructed values
    ϕᵢ₊½ = (ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i]) + ϕᴿᵢ₊½(ϕ[ᵢ₊₁], ∂ϕ[ᵢ₊₁], ∂²ϕ[ᵢ₊₁])) / 2
    ϕᵢ₋½ = (ϕᴸᵢ₊½(ϕ[ᵢ₋₁], ∂ϕ[ᵢ₋₁], ∂²ϕ[ᵢ₋₁]) + ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])) / 2

    _ϕ_ξ = ϕᵢ₊½ - ϕᵢ₋½
    _ϕ_ξ = _ϕ_ξ * (abs(_ϕ_ξ) >= ϵ)
    ϕ_ξ[i] = _ϕ_ξ
  end

  # if !inner_domain_only
  index_offset = 0

  # select first index along given boundary axis
  lo_domain = lower_boundary_indices(domain, axis, index_offset)

  # the lo-side derivative ∂ϕ/∂ξ can only use ϕᴿᵢ₋½ instead of ϕᵢ₋½
  for i in lo_domain
    ᵢ₊₁ = shift(i, axis, +1)

    ϕᵢ₊½ = (ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i]) + ϕᴿᵢ₊½(ϕ[ᵢ₊₁], ∂ϕ[ᵢ₊₁], ∂²ϕ[ᵢ₊₁])) / 2
    ϕᵢ₋½ = ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])

    _ϕ_ξ = ϕᵢ₊½ - ϕᵢ₋½
    _ϕ_ξ = _ϕ_ξ * (abs(_ϕ_ξ) >= ϵ)
    ϕ_ξ[i] = _ϕ_ξ
  end

  # select last index along given boundary axis
  hi_domain = upper_boundary_indices(domain, axis, index_offset)

  # the hi-side derivative ∂ϕ/∂ξ can only use ϕᴸᵢ₊½ instead of ϕᵢ₊½
  for i in hi_domain
    ᵢ₋₁ = shift(i, axis, -1)

    ϕᵢ₊½ = ϕᴸᵢ₊½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])
    ϕᵢ₋½ = (ϕᴸᵢ₊½(ϕ[ᵢ₋₁], ∂ϕ[ᵢ₋₁], ∂²ϕ[ᵢ₋₁]) + ϕᴿᵢ₋½(ϕ[i], ∂ϕ[i], ∂²ϕ[i])) / 2

    _ϕ_ξ = ϕᵢ₊½ - ϕᵢ₋½
    _ϕ_ξ = _ϕ_ξ * (abs(_ϕ_ξ) >= ϵ)
    ϕ_ξ[i] = _ϕ_ξ
  end
  # end

  return nothing
end

function first_deriv!(∂ϕ::AbstractArray, ϕ::AbstractArray, axis::Int, domain, nhalo=3)
  # are we computing the derivatives on the inner portion of the domain only? if so
  # we don't need to use one-sided stencils along the edges
  # if size(domain) == size(ϕ)
  # inner_domain_only = false
  # domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, axis, -nhalo)
  # else
  #   inner_domain_only = true
  #   inner_domain = domain
  # end

  for i in inner_domain
    ϕᵢ₋₃ = ϕ[shift(i, axis, -3)]
    ϕᵢ₋₂ = ϕ[shift(i, axis, -2)]
    ϕᵢ₋₁ = ϕ[shift(i, axis, -1)]
    ϕᵢ₊₁ = ϕ[shift(i, axis, +1)]
    ϕᵢ₊₂ = ϕ[shift(i, axis, +2)]
    ϕᵢ₊₃ = ϕ[shift(i, axis, +3)]

    ∂ϕ[i] = ∂(ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃)
  end

  # if !inner_domain_only
  index_offset = 0

  # select first index along given boundary axis
  lo_domain = lower_boundary_indices(domain, axis, index_offset)

  # and then iterate along that domain to calculate the derivatives
  # of the first three cells using one-sided and mixed-offset stencils
  for i in lo_domain
    ᵢ₊₁ = shift(i, axis, +1)
    ᵢ₊₂ = shift(i, axis, +2)
    ᵢ₊₃ = shift(i, axis, +3)
    ᵢ₊₄ = shift(i, axis, +4)
    ᵢ₊₅ = shift(i, axis, +5)

    ∂ϕ₁, ∂ϕ₂, ∂ϕ₃ = ∂_loedge_6th_order(ϕ[i], ϕ[ᵢ₊₁], ϕ[ᵢ₊₂], ϕ[ᵢ₊₃], ϕ[ᵢ₊₄], ϕ[ᵢ₊₅])

    ∂ϕ[i] = ∂ϕ₁
    ∂ϕ[ᵢ₊₁] = ∂ϕ₂
    ∂ϕ[ᵢ₊₂] = ∂ϕ₃
  end

  # select last index along given boundary axis
  hi_domain = upper_boundary_indices(domain, axis, index_offset)

  # and then iterate along that domain to calculate the derivatives
  # of the last three cells using one-sided and mixed-offset stencils
  for i in hi_domain
    ᵢ₋₁ = shift(i, axis, -1)
    ᵢ₋₂ = shift(i, axis, -2)
    ᵢ₋₃ = shift(i, axis, -3)
    ᵢ₋₄ = shift(i, axis, -4)
    ᵢ₋₅ = shift(i, axis, -5)

    ∂ϕ₃, ∂ϕ₂, ∂ϕ₁ = ∂_hiedge_6th_order(ϕ[ᵢ₋₅], ϕ[ᵢ₋₄], ϕ[ᵢ₋₃], ϕ[ᵢ₋₂], ϕ[ᵢ₋₁], ϕ[i])

    ∂ϕ[i] = ∂ϕ₁
    ∂ϕ[ᵢ₋₁] = ∂ϕ₂
    ∂ϕ[ᵢ₋₂] = ∂ϕ₃
  end
  # end
end

@kernel inbounds = true function _first_deriv_inner_kernel(∂ϕ, ϕ, axis, domain)
  idx = @index(Global, Linear)
  i = domain[idx]
  ϕᵢ₋₃ = ϕ[shift(i, axis, -3)]
  ϕᵢ₋₂ = ϕ[shift(i, axis, -2)]
  ϕᵢ₋₁ = ϕ[shift(i, axis, -1)]
  ϕᵢ₊₁ = ϕ[shift(i, axis, +1)]
  ϕᵢ₊₂ = ϕ[shift(i, axis, +2)]
  ϕᵢ₊₃ = ϕ[shift(i, axis, +3)]

  ∂ϕ[i] = ∂(ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃)
end

function second_deriv!(
  ∂²ϕ::AbstractArray, ∂ϕ::AbstractArray, ϕ::AbstractArray, axis::Int, domain, nhalo=3
)

  # are we computing the derivatives on the inner portion of the domain only? if so
  # we don't need to use one-sided stencils along the edges
  # if size(domain) == size(ϕ)
  #   inner_domain_only = false
  # domain = CartesianIndices(ϕ)
  inner_domain = expand(domain, axis, -nhalo)
  # else
  #   inner_domain_only = true
  #   inner_domain = domain
  # end

  for i in inner_domain
    ᵢ₋₁ = shift(i, axis, -1)
    ᵢ₊₁ = shift(i, axis, +1)

    ∂²ϕ[i] = ∂²(∂ϕ[ᵢ₊₁], ∂ϕ[ᵢ₋₁], ϕ[ᵢ₋₁], ϕ[i], ϕ[ᵢ₊₁])
  end

  # The edges require one-sided or mixed stencils. For the 2nd derivative, 
  # just call take the deriv of the deriv
  # if !inner_domain_only
  index_offset = 0

  # select first index along given boundary axis
  lo_domain = lower_boundary_indices(domain, axis, index_offset)

  # and then iterate along that domain to calculate the derivatives
  # of the first three cells using one-sided and mixed-offset stencils
  for i in lo_domain
    ᵢ₊₁ = shift(i, axis, +1)
    ᵢ₊₂ = shift(i, axis, +2)
    ᵢ₊₃ = shift(i, axis, +3)
    ᵢ₊₄ = shift(i, axis, +4)
    ᵢ₊₅ = shift(i, axis, +5)

    ∂²ϕ₁, ∂²ϕ₂, ∂²ϕ₃ = ∂_loedge_6th_order(
      ∂ϕ[i], ∂ϕ[ᵢ₊₁], ∂ϕ[ᵢ₊₂], ∂ϕ[ᵢ₊₃], ∂ϕ[ᵢ₊₄], ∂ϕ[ᵢ₊₅]
    )

    ∂²ϕ[i] = ∂²ϕ₁
    ∂²ϕ[ᵢ₊₁] = ∂²ϕ₂
    ∂²ϕ[ᵢ₊₂] = ∂²ϕ₃
  end

  # select last index along given boundary axis
  hi_domain = upper_boundary_indices(domain, axis, index_offset)

  # and then iterate along that domain to calculate the derivatives
  # of the last three cells using one-sided and mixed-offset stencils
  for i in hi_domain
    ᵢ₋₁ = shift(i, axis, -1)
    ᵢ₋₂ = shift(i, axis, -2)
    ᵢ₋₃ = shift(i, axis, -3)
    ᵢ₋₄ = shift(i, axis, -4)
    ᵢ₋₅ = shift(i, axis, -5)

    ∂²ϕ₃, ∂²ϕ₂, ∂²ϕ₁ = ∂_hiedge_6th_order(
      ∂ϕ[ᵢ₋₅], ∂ϕ[ᵢ₋₄], ∂ϕ[ᵢ₋₃], ∂ϕ[ᵢ₋₂], ∂ϕ[ᵢ₋₁], ∂ϕ[i]
    )

    ∂²ϕ[i] = ∂²ϕ₁
    ∂²ϕ[ᵢ₋₁] = ∂²ϕ₂
    ∂²ϕ[ᵢ₋₂] = ∂²ϕ₃
  end
  # end
end

@kernel inbounds = true function _second_deriv_kernel!(∂²ϕ, ∂ϕ, ϕ, axis, domain)
  idx = @index(Global, Linear)
end