
struct MonotoneExplicitGradientScheme{N,C,DS} <: DiscretizationScheme
  cache::C
  nhalo::Int
  derivative_scheme::DS
  use_symmetric_conservative_metric_scheme::Bool
end

function Adapt.adapt_structure(
  to, scheme::MonotoneExplicitGradientScheme{N,C,DS}
) where {N,C<:Nothing,DS}

  #
  return MonotoneExplicitGradientScheme{N,Nothing,DS}(
    nothing,
    scheme.nhalo,
    scheme.derivative_scheme,
    scheme.use_symmetric_conservative_metric_scheme,
  )
end

function Adapt.adapt_structure(
  to, scheme::MonotoneExplicitGradientScheme{N,C,DS}
) where {N,C,DS}

  #
  cache = (;
    ∂²ϕ=Adapt.adapt_structure(to, scheme.cache.∂²ϕ),
    ∂ϕ=Adapt.adapt_structure(to, scheme.cache.∂ϕ),
    outer_deriv_1=Adapt.adapt_structure(to, scheme.cache.outer_deriv_1),
    outer_deriv_2=Adapt.adapt_structure(to, scheme.cache.outer_deriv_2),
    inner_deriv_1=Adapt.adapt_structure(to, scheme.cache.inner_deriv_1),
    inner_deriv_2=Adapt.adapt_structure(to, scheme.cache.inner_deriv_2),
  )
  return MonotoneExplicitGradientScheme{N,typeof(cache),DS}(
    cache,
    scheme.nhalo,
    scheme.derivative_scheme,
    scheme.use_symmetric_conservative_metric_scheme,
  )
end

# The MEG6 scheme requires a halo of 5 cells in all dimensions

"""
    MonotoneExplicitGradientScheme(order::Int; use_cache=true, celldims=nothing, backend=CPU(), T=Float64, use_symmetric_conservative_metric_scheme=false)

Create a monotone explicit gradient (MEG) discretization of the requested `order` (2, 4, or 6).

When `use_cache` is `true`, temporary derivative buffers sized by `celldims` are preallocated on the chosen `backend` using element type `T`. The `use_symmetric_conservative_metric_scheme` flag enables symmetric conservative metric updates when metrics are refreshed.
"""
function MonotoneExplicitGradientScheme(
  order::Int;
  use_cache=true,
  celldims=nothing,
  backend=CPU(),
  T=Float64,
  use_symmetric_conservative_metric_scheme=false,
)
  if order == 2
    nhalo = nhalo_lookup[:MEG2] # 2
    deriv_scheme = SecondOrder()
  elseif order == 4
    nhalo = nhalo_lookup[:MEG4] # 3
    deriv_scheme = FourthOrder()
  elseif order == 6
    nhalo = nhalo_lookup[:MEG6] # 5
    deriv_scheme = SixthOrder()
  else
    error("Unsupported derivative order")
  end

  if use_cache && isnothing(celldims)
    error(
      "When use_cache=true (default), celldims must be an NTuple providing the dimensions of the cell-based domain",
    )
  end

  if !isnothing(celldims)
    if !all((celldims .- 2nhalo) .> 5)
      @warn "The domain dimensions $(celldims) are too small to use a 6th order scheme"
    end
  end

  if use_cache
    cache = (;
      ∂²ϕ=KernelAbstractions.zeros(backend, T, celldims),
      ∂ϕ=KernelAbstractions.zeros(backend, T, celldims),
      outer_deriv_1=KernelAbstractions.zeros(backend, T, celldims),
      outer_deriv_2=KernelAbstractions.zeros(backend, T, celldims),
      inner_deriv_1=KernelAbstractions.zeros(backend, T, celldims),
      inner_deriv_2=KernelAbstractions.zeros(backend, T, celldims),
    )
  else
    cache = nothing
  end

  return MonotoneExplicitGradientScheme{order,typeof(cache),typeof(deriv_scheme)}(
    cache, nhalo, deriv_scheme, use_symmetric_conservative_metric_scheme
  )
end

"""
    ϕᴸᵢ₊½(ϕᵢ, ∂ϕᵢ, ∂²ϕᵢ)

Get the reconstructed left edge value at i+1/2 using the cell-center value `ϕᵢ`, derivative `∂ϕᵢ`, and second derivative `∂²ϕᵢ`
"""
@inline ϕᴸᵢ₊½(ϕᵢ, ∂ϕᵢ, ∂²ϕᵢ) = ϕᵢ + ((1 / 2) * ∂ϕᵢ + (1 / 12) * ∂²ϕᵢ)

"""
    ϕᴿᵢ₊½(ϕᵢ₊₁, ∂ϕᵢ₊₁, ∂²ϕᵢ₊₁)

Get the reconstructed right edge value at i+1/2 using the cell-center value `ϕᵢ`, derivative `∂ϕᵢ`, and second derivative `∂²ϕᵢ`
"""
@inline ϕᴿᵢ₊½(ϕᵢ₊₁, ∂ϕᵢ₊₁, ∂²ϕᵢ₊₁) = ϕᵢ₊₁ - ((1 / 2) * ∂ϕᵢ₊₁ + (1 / 12) * ∂²ϕᵢ₊₁)

"""
    ϕᴿᵢ₋½(ϕᵢ, ∂ϕᵢ, ∂²ϕᵢ)

Get the reconstructed right edge value at i-1/2 using the cell-center value `ϕᵢ`, derivative `∂ϕᵢ`, and second derivative `∂²ϕᵢ`
"""
@inline ϕᴿᵢ₋½(ϕᵢ, ∂ϕᵢ, ∂²ϕᵢ) = ϕᵢ - ((1 / 2) * ∂ϕᵢ + (1 / 12) * ∂²ϕᵢ)
