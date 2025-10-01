
function central_first_derivative!(order, ∂W, W, axis::Int, domain, backend)
  central_first_derivative_inner_domain_kernel!(order, ∂W, W, axis, domain, backend)
end

function central_second_derivative!(order, ∂²W, ∂W, W, axis::Int, domain, backend)
  second_derivative_inner_domain_kernel!(∂²W, ∂W, W, axis, domain, backend)
end

# ------------------------
# Central First Derivative
# ------------------------

## CPU Kernels

function central_first_derivative_inner_domain_kernel!(
  ::SecondOrder, ∂W, W, axis, domain, ::CPU
)

  #
  @batch for I in domain
    ᵢ₊₁ = shift(I, axis, +1)
    ᵢ₋₁ = shift(I, axis, -1)

    @inline ∂W[I] = central_derivative(W[ᵢ₋₁], W[ᵢ₊₁])
  end
end

function central_first_derivative_inner_domain_kernel!(
  ::FourthOrder, ∂W, W, axis, domain, ::CPU
)

  #
  @batch for I in domain
    ᵢ₊₁ = shift(I, axis, +1)
    ᵢ₊₂ = shift(I, axis, +2)
    ᵢ₋₁ = shift(I, axis, -1)
    ᵢ₋₂ = shift(I, axis, -2)

    @inline ∂W[I] = central_derivative(W[ᵢ₋₂], W[ᵢ₋₁], W[ᵢ₊₁], W[ᵢ₊₂])
  end
end

function central_first_derivative_inner_domain_kernel!(
  ::SixthOrder, ∂W, W, axis, domain, ::CPU
)

  #
  @batch for I in domain
    ᵢ₊₁ = shift(I, axis, +1)
    ᵢ₊₂ = shift(I, axis, +2)
    ᵢ₊₃ = shift(I, axis, +3)
    ᵢ₋₁ = shift(I, axis, -1)
    ᵢ₋₂ = shift(I, axis, -2)
    ᵢ₋₃ = shift(I, axis, -3)

    @inline ∂W[I] = central_derivative(W[ᵢ₋₃], W[ᵢ₋₂], W[ᵢ₋₁], W[ᵢ₊₁], W[ᵢ₊₂], W[ᵢ₊₃])
  end
end

function central_first_derivative_inner_domain_kernel!(
  ::EighthOrder, ∂W, W, axis, domain, ::CPU
)

  #
  @batch for I in domain
    ᵢ₊₁ = shift(I, axis, +1)
    ᵢ₊₂ = shift(I, axis, +2)
    ᵢ₊₃ = shift(I, axis, +3)
    ᵢ₊₄ = shift(I, axis, +4)
    ᵢ₋₁ = shift(I, axis, -1)
    ᵢ₋₂ = shift(I, axis, -2)
    ᵢ₋₃ = shift(I, axis, -3)
    ᵢ₋₄ = shift(I, axis, -4)

    @inline ∂W[I] = central_derivative(
      W[ᵢ₋₄], W[ᵢ₋₃], W[ᵢ₋₂], W[ᵢ₋₁], W[ᵢ₊₁], W[ᵢ₊₂], W[ᵢ₊₃], W[ᵢ₊₄]
    )
  end
end

## GPU Kernels

function central_first_derivative_inner_domain_kernel!(
  ::SecondOrder, ∂W, W, axis, domain, backend::GPU
)

  #
  @kernel inbounds = true function _kernel(∂ϕ, ϕ, iaxis, ∂domain)
    idx = @index(Global, Linear)
    I = ∂domain[idx]

    ᵢ₊₁ = shift(I, iaxis, +1)
    ᵢ₋₁ = shift(I, iaxis, -1)

    @inline ∂ϕ[I] = central_derivative(ϕ[ᵢ₋₁], ϕ[ᵢ₊₁])
  end

  _kernel(backend)(∂W, W, iaxis, ∂domain; ndrange=size(∂domain))

  return nothing
end

function central_first_derivative_inner_domain_kernel!(
  ::FourthOrder, ∂W, W, axis, domain, backend::GPU
)

  #
  @kernel inbounds = true function _kernel(∂ϕ, ϕ, iaxis, ∂domain)
    idx = @index(Global, Linear)
    I = ∂domain[idx]

    ᵢ₊₁ = shift(I, iaxis, +1)
    ᵢ₊₂ = shift(I, iaxis, +2)
    ᵢ₋₁ = shift(I, iaxis, -1)
    ᵢ₋₂ = shift(I, iaxis, -2)

    @inline ∂ϕ[I] = central_derivative(ϕ[ᵢ₋₂], ϕ[ᵢ₋₁], ϕ[ᵢ₊₁], ϕ[ᵢ₊₂])
  end

  _kernel(backend)(∂W, W, iaxis, ∂domain; ndrange=size(∂domain))

  return nothing
end

function central_first_derivative_inner_domain_kernel!(
  ::SixthOrder, ∂W, W, axis, domain, backend::GPU
)

  #
  @kernel inbounds = true function _kernel(∂ϕ, ϕ, iaxis, ∂domain)
    idx = @index(Global, Linear)
    I = ∂domain[idx]

    ᵢ₊₁ = shift(I, iaxis, +1)
    ᵢ₊₂ = shift(I, iaxis, +2)
    ᵢ₊₃ = shift(I, iaxis, +3)
    ᵢ₋₁ = shift(I, iaxis, -1)
    ᵢ₋₂ = shift(I, iaxis, -2)
    ᵢ₋₃ = shift(I, iaxis, -3)

    @inline ∂ϕ[I] = central_derivative(ϕ[ᵢ₋₃], ϕ[ᵢ₋₂], ϕ[ᵢ₋₁], ϕ[ᵢ₊₁], ϕ[ᵢ₊₂], ϕ[ᵢ₊₃])
  end

  _kernel(backend)(∂W, W, iaxis, ∂domain; ndrange=size(∂domain))

  return nothing
end

function central_first_derivative_inner_domain_kernel!(
  ::EighthOrder, ∂W, W, axis, domain, backend::GPU
)

  #
  @kernel inbounds = true function _kernel(∂ϕ, ϕ, iaxis, ∂domain)
    idx = @index(Global, Linear)
    I = ∂domain[idx]

    ᵢ₊₁ = shift(I, iaxis, +1)
    ᵢ₊₂ = shift(I, iaxis, +2)
    ᵢ₊₃ = shift(I, iaxis, +3)
    ᵢ₊₄ = shift(I, iaxis, +4)
    ᵢ₋₁ = shift(I, iaxis, -1)
    ᵢ₋₂ = shift(I, iaxis, -2)
    ᵢ₋₃ = shift(I, iaxis, -3)
    ᵢ₋₄ = shift(I, iaxis, -4)

    @inline ∂ϕ[I] = central_derivative(
      ϕ[ᵢ₋₄], ϕ[ᵢ₋₃], ϕ[ᵢ₋₂], ϕ[ᵢ₋₁], ϕ[ᵢ₊₁], ϕ[ᵢ₊₂], ϕ[ᵢ₊₃], ϕ[ᵢ₊₄]
    )
  end

  _kernel(backend)(∂W, W, iaxis, ∂domain; ndrange=size(∂domain))

  return nothing
end

# ---------------------------------
# Mixed First Derivative for edges
# ---------------------------------

## CPU Kernels

function mixed_first_derivative_lo_edge_kernel!(
  ::SixthOrder, ∂W, W, axis, domain, backend::CPU
)
  index_offset = 0
  lo_domain = lower_boundary_indices(domain, axis, index_offset)

  @batch for i in lo_domain
    ᵢ₊₁ = shift(i, axis, +1)
    ᵢ₊₂ = shift(i, axis, +2)
    ᵢ₊₃ = shift(i, axis, +3)
    ᵢ₊₄ = shift(i, axis, +4)
    ᵢ₊₅ = shift(i, axis, +5)

    ∂W₁, ∂W₂, ∂W₃ = mixed_derivatives_loedge(W[i], W[ᵢ₊₁], W[ᵢ₊₂], W[ᵢ₊₃], W[ᵢ₊₄], W[ᵢ₊₅])

    ∂W[i] = ∂W₁
    ∂W[ᵢ₊₁] = ∂W₂
    ∂W[ᵢ₊₂] = ∂W₃
  end
end

function mixed_first_derivative_hi_edge_kernel!(
  ::SixthOrder, ∂W, W, axis, domain, backend::CPU
)
  index_offset = 0
  hi_domain = upper_boundary_indices(domain, axis, index_offset)

  @batch for i in hi_domain
    ᵢ₋₁ = shift(i, axis, -1)
    ᵢ₋₂ = shift(i, axis, -2)
    ᵢ₋₃ = shift(i, axis, -3)
    ᵢ₋₄ = shift(i, axis, -4)
    ᵢ₋₅ = shift(i, axis, -5)

    ∂W₃, ∂W₂, ∂W₁ = mixed_derivatives_hiedge(W[ᵢ₋₅], W[ᵢ₋₄], W[ᵢ₋₃], W[ᵢ₋₂], W[ᵢ₋₁], W[i])

    ∂W[i] = ∂W₁
    ∂W[ᵢ₋₁] = ∂W₂
    ∂W[ᵢ₋₂] = ∂W₃
  end
end

## GPU Kernels
function mixed_first_derivative_lo_edge_kernel!(
  ::SixthOrder, ∂W, W, ∂axis, domain, backend::GPU
)
  @kernel inbounds = true function _first_deriv_lo_kernel!(ϕ, ∂ϕ, axis, ∂domain)
    idx = @index(Global, Linear)
    i = ∂domain[idx]

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

  # select first index along given boundary axis
  lo_domain = lower_boundary_indices(domain, ∂axis, index_offset)

  # and then iterate along that domain to calculate the derivatives
  # of the first three cells using one-sided and mixed-offset stencils
  _first_deriv_lo_kernel!(backend)(W, ∂W, ∂axis, lo_domain; ndrange=size(lo_domain))
end

function mixed_first_derivative_hi_edge_kernel!(
  ::SixthOrder, ∂W, W, ∂axis, domain, backend::GPU
)
  @kernel inbounds = true function _first_deriv_hi_kernel!(ϕ, ∂ϕ, axis, ∂domain)
    idx = @index(Global, Linear)
    i = ∂domain[idx]

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

  # select last index along given boundary axis
  hi_domain = upper_boundary_indices(domain, ∂axis, index_offset)

  # and then iterate along that domain to calculate the derivatives
  # of the last three cells using one-sided and mixed-offset stencils
  _first_deriv_hi_kernel!(backend)(ϕ, ∂ϕ, axis, hi_domain; ndrange=size(hi_domain))
end

# ------------------------
# Second Derivative
# ------------------------
function central_second_derivative_inner_domain_kernel!(∂²W, ∂W, W, axis, domain, ::CPU)
  @batch for I in domain
    ᵢ₊₁ = shift(I, axis, +1)
    ᵢ₋₁ = shift(I, axis, -1)

    @inline ∂²W[I] = second_derivative(W[ᵢ₋₁], W[I], W[ᵢ₊₁], ∂W[ᵢ₋₁], ∂W[ᵢ₊₁])
  end
end

function central_second_derivative_inner_domain_kernel!(
  ∂²W, ∂W, W, axis, domain, backend::GPU
)

  #
  @kernel inbounds = true function _kernel(∂²W, ∂ϕ, ϕ, iaxis, ∂domain)
    idx = @index(Global, Linear)
    I = ∂domain[idx]
    ᵢ₊₁ = shift(I, iaxis, +1)
    ᵢ₋₁ = shift(I, iaxis, -1)

    @inline ∂²W[I] = second_derivative(ϕ[ᵢ₋₁], ϕ[I], ϕ[ᵢ₊₁], ∂ϕ[ᵢ₋₁], ∂ϕ[ᵢ₊₁])
  end

  _kernel(backend)(∂²W, ∂W, W, axis, domain; ndrange=size(domain))

  return nothing
end