module DiscretizationSchemes

using StaticArrays
using KernelAbstractions
using Polyester
using MappedArrays
using CartesianDomains
using Adapt
using LinearAlgebra

export MonotoneExplicitGradientScheme

export SecondOrder, FourthOrder, SixthOrder, EighthOrder
export compute_first_derivatives!, compute_second_derivatives!
export cell_center_derivatives!, interpolate_to_edge!

abstract type DiscretizationScheme end
abstract type DerivativeScheme end

"""
    SecondOrder()

Marker type indicating a second-order centered derivative stencil.
"""
struct SecondOrder <: DerivativeScheme end
"""
    FourthOrder()

Marker type indicating a fourth-order centered derivative stencil.
"""
struct FourthOrder <: DerivativeScheme end
"""
    SixthOrder()

Marker type indicating a sixth-order centered derivative stencil.
"""
struct SixthOrder <: DerivativeScheme end
"""
    EighthOrder()

Marker type indicating an eighth-order centered derivative stencil.
"""
struct EighthOrder <: DerivativeScheme end

include("derivative_stencils.jl")
include("derivative_kernels.jl")

include("meg/MonotoneExplicitGradients.jl")
include("meg/cell_center_derivs.jl")
include("meg/interp_to_edges.jl")
# first_derivative_domain(scheme, domain) = expand(domain, -scheme.nhalo ÷ 2)

const nhalo_lookup = Dict(
  :MEG2 => 2,
  :MEG4 => 3,
  :MEG6 => 5,
  :MEG2_SYMMETRIC => 2,
  :MEG4_SYMMETRIC => 3,
  :MEG6_SYMMETRIC => 5,
)

#! format: off
first_derivative_domain(::MonotoneExplicitGradientScheme{2}, domain) = expand(domain, -1)
first_derivative_domain(::MonotoneExplicitGradientScheme{4}, domain) = expand(domain, -2)
first_derivative_domain(::MonotoneExplicitGradientScheme{6}, domain) = expand(domain, -3)
first_derivative_domain(::MonotoneExplicitGradientScheme{2}, domain, axis) = expand(domain, axis, -1)
first_derivative_domain(::MonotoneExplicitGradientScheme{4}, domain, axis) = expand(domain, axis, -2)
first_derivative_domain(::MonotoneExplicitGradientScheme{6}, domain, axis) = expand(domain, axis, -3)

# second_derivative_domain(scheme, domain)
second_derivative_domain(::MonotoneExplicitGradientScheme{2}, domain) = expand(domain, -2)
second_derivative_domain(::MonotoneExplicitGradientScheme{4}, domain) = expand(domain, -3)
second_derivative_domain(::MonotoneExplicitGradientScheme{6}, domain) = expand(domain, -4)
second_derivative_domain(::MonotoneExplicitGradientScheme{2}, domain, axis) = expand(domain, axis, -2)
second_derivative_domain(::MonotoneExplicitGradientScheme{4}, domain, axis) = expand(domain, axis, -3)
second_derivative_domain(::MonotoneExplicitGradientScheme{6}, domain, axis) = expand(domain, axis, -4)

inner_domain_with_one_sided_edges(::MonotoneExplicitGradientScheme{2}, ∂domain, axis) = expand(∂domain, axis, -1)
inner_domain_with_one_sided_edges(::MonotoneExplicitGradientScheme{4}, ∂domain, axis) = expand(∂domain, axis, -2)
inner_domain_with_one_sided_edges(::MonotoneExplicitGradientScheme{6}, ∂domain, axis) = expand(∂domain, axis, -3)

#! format: on
"""
    compute_first_derivatives!(scheme::DiscretizationScheme, ∂ϕ, ϕ, axis, derivative_domain, backend, use_one_sided_on_edges)

Compute the 1st derivative `∂ϕ` of `ϕ` () along the given `axis`. The order of the scheme determines how many padded halo cells around the edges are to be used. When `use_one_sided_on_edges`
is `true`, one-sided derivatives will be used along the edges. If `false`, central derivatives will be everywhere, so if the data
within the halo region is non-sensical, it will produce bad data.
"""
function compute_first_derivatives!(
  scheme::DiscretizationScheme,
  ∂ϕ,
  ϕ::AbstractArray{T,N},
  axis,
  domain,
  backend,
  use_one_sided_on_edges::Bool,
  floor=eps(T),
) where {T,N}
  if use_one_sided_on_edges
    # shrink the derivative domain so the central deriv form can be used
    inner_domain = inner_domain_with_one_sided_edges(scheme, domain, axis)
  else
    # otherwise, use the central deriv form everywhere
    inner_domain = domain
  end

  central_first_derivative_inner_domain_kernel!(
    scheme.derivative_scheme, ∂ϕ, ϕ, axis, inner_domain, backend
  )

  if use_one_sided_on_edges
    mixed_first_derivative_lo_edge_kernel!(
      scheme.derivative_scheme, ∂ϕ, ϕ, axis, domain, backend
    )
    mixed_first_derivative_hi_edge_kernel!(
      scheme.derivative_scheme, ∂ϕ, ϕ, axis, domain, backend
    )
  end
end

"""
    compute_second_derivatives!(scheme::DiscretizationScheme, ∂²ϕ, ∂ϕ, ϕ, axis, derivative_domain, backend, use_one_sided_on_edges)

Compute the 2nd derivative `∂²ϕ` of `ϕ` () along the given `axis`. The order of the scheme determines how many padded halo cells around the edges are to be used. When `use_one_sided_on_edges`
is `true`, one-sided derivatives will be used along the edges. If `false`, central derivatives will be everywhere, so if the data
within the halo region is non-sensical, it will produce bad data.
"""
function compute_second_derivatives!(
  scheme::DiscretizationScheme,
  ∂²ϕ,
  ∂ϕ,
  ϕ::AbstractArray{T,N},
  axis,
  derivative_domain,
  backend,
  use_one_sided_on_edges::Bool,
  floor=eps(T),
) where {T,N}
  if use_one_sided_on_edges
    # shrink the domain so the central deriv form can be used
    # both ϕ and ∂ϕ in this form are _not_ populated in the halo regions
    central_domain = inner_domain_with_one_sided_edges(scheme, derivative_domain, axis)
  else
    # both ϕ and ∂ϕ need to be populated in the halo regions for this to work properly
    central_domain = derivative_domain
  end

  central_second_derivative_inner_domain_kernel!(∂²ϕ, ∂ϕ, ϕ, axis, central_domain, backend)

  if use_one_sided_on_edges
    mixed_first_derivative_lo_edge_kernel!(
      scheme.derivative_scheme, ∂²ϕ, ∂ϕ, axis, derivative_domain, backend
    )
    mixed_first_derivative_hi_edge_kernel!(
      scheme.derivative_scheme, ∂²ϕ, ∂ϕ, axis, derivative_domain, backend
    )
  end
end

end
