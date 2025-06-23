module MontoneExplicitGradientSchemes

using StaticArrays
using KernelAbstractions
using MappedArrays
using CartesianDomains
using MappedArrays

using ..DiscretizationSchemes
using ..DiscretizationSchemes: DiscretizationScheme

export MontoneExplicitGradientScheme6thOrder
export reconstruct_interface,
  interpolate_to_edge,
  interface_derivative,
  interface_derivative_nodamping,
  cell_center_derivative

export cell_center_derivatives!, interpolate_to_edge!

struct MontoneExplicitGradientScheme6thOrder{C,B} <: DiscretizationScheme
  cache::C
  backend::B
  nhalo::Int
  use_symmetric_conservative_metric_scheme::Bool
end

# The MEG6 scheme requires a halo of 5 cells in all dimensions

function MontoneExplicitGradientScheme6thOrder(;
  use_cache=true,
  celldims=nothing,
  backend=CPU(),
  T=Float64,
  use_symmetric_conservative_metric_scheme=false,
)
  nhalo = 5
  if use_cache && isnothing(celldims)
    error(
      "When use_cache=true (default), celldims must be an NTuple providing the dimensions of the cell-based domain",
    )
  end

  if !all((celldims .- 2nhalo) .> 5)
    @warn "The domain dimensions $(celldims) are too small to use a 6th order scheme"
  end

  if use_cache
    cache = (;
      ∂²ϕ=KernelAbstractions.zeros(backend, T, celldims),
      ∂ϕ=KernelAbstractions.zeros(backend, T, celldims),
      outer_deriv_1=KernelAbstractions.zeros(backend, T, celldims),
      outer_deriv_2=KernelAbstractions.zeros(backend, T, celldims),
    )
  else
    cache = nothing
  end

  return MontoneExplicitGradientScheme6thOrder(
    cache, backend, nhalo, use_symmetric_conservative_metric_scheme
  )
end

include("cell_center_derivs.jl")
include("interp_to_edges.jl")

reconstruct_interface(::MontoneExplicitGradientScheme6thOrder, ϕ...) = meg6_ϕᴸᴿᵢ₊½(ϕ...)

interpolate_to_edge(::MontoneExplicitGradientScheme6thOrder, ϕ...) = meg6_ϕᵢ₊½(ϕ...)

function interface_derivative(::MontoneExplicitGradientScheme6thOrder, ϕ...)
  inner_interface_deriv_6th_order_alpha_damping(ϕ...)
end

function interface_derivative_nodamping(::MontoneExplicitGradientScheme6thOrder, ϕ...)
  inner_interface_deriv_6th_order(ϕ...)
end

"""
    cell_center_derivative(::MontoneExplicitGradientScheme6thOrder, ϕᵢ₋₅, ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅)

TBW
"""
function DiscretizationSchemes.cell_center_derivative(
  ::MontoneExplicitGradientScheme6thOrder,
  ϕᵢ₋₅,
  ϕᵢ₋₄,
  ϕᵢ₋₃,
  ϕᵢ₋₂,
  ϕᵢ₋₁,
  ϕᵢ₊₁,
  ϕᵢ₊₂,
  ϕᵢ₊₃,
  ϕᵢ₊₄,
  ϕᵢ₊₅,
)
  inner_central_deriv_6th_order(ϕᵢ₋₅, ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅)
end

# 2nd order gradient
function ∂(ϕᵢ₋₁, ϕᵢ₊₁)
  ∂ϕᵢ = (ϕᵢ₊₁ - ϕᵢ₋₁) / 2
  return ∂ϕᵢ
end

# 4th order gradient
function ∂(ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂)
  ∂ϕᵢ = (
    (2 / 3) * (ϕᵢ₊₁ - ϕᵢ₋₁) + # 
    (1 / 12) * (ϕᵢ₋₂ - ϕᵢ₊₂)  #
  )
  return ∂ϕᵢ
end

# 6th order gradient
function ∂(ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃)
  ∂ϕᵢ = (
    (3 / 4) * (ϕᵢ₊₁ - ϕᵢ₋₁) -  #
    (3 / 20) * (ϕᵢ₊₂ - ϕᵢ₋₂) + # 
    (1 / 60) * (ϕᵢ₊₃ - ϕᵢ₋₃)   #
  )

  return ∂ϕᵢ
end

function ∂_loedge_6th_order(ϕ₁, ϕ₂, ϕ₃, ϕ₄, ϕ₅, ϕ₆)
  # FiniteDifferenceMethod(0:5, 1)
  ∂ϕ₁ = (
    (-137 / 60) * ϕ₁ + #
    5ϕ₂ -            #
    5ϕ₃ +            #
    (10 / 3) * ϕ₄ -  #
    (5 / 4) * ϕ₅ +   #
    (1 / 5) * ϕ₆     #
  )

  # FiniteDifferenceMethod(-1:4, 1)
  ∂ϕ₂ = (
    (-1 / 5) * ϕ₁ -    #
    (13 / 12) * ϕ₂ + #
    2ϕ₃ -            #
    ϕ₄ +             #
    (1 / 3) * ϕ₅ -   #
    (1 / 20) * ϕ₆    #
  )

  # FiniteDifferenceMethod(-2:3, 1)
  ∂ϕ₃ = (
    (1 / 20) * ϕ₁ -   #
    (1 / 2) * ϕ₂ -  #
    (1 / 3) * ϕ₃ +  #
    ϕ₄ -            #
    (1 / 4) * ϕ₅ +  #
    (1 / 30) * ϕ₆   #
  )

  return (∂ϕ₁, ∂ϕ₂, ∂ϕ₃)
end

function ∂_hiedge_6th_order(ϕ₁, ϕ₂, ϕ₃, ϕ₄, ϕ₅, ϕ₆)
  # FiniteDifferenceMethod(-5:0, 1)
  ∂ϕ₆ = (
    (-1 / 5) * ϕ₁ +   #
    (5 / 4) * ϕ₂ -    #
    (10 / 3) * ϕ₃ +   #
    5ϕ₄ -             #
    5ϕ₅ +             #
    (137//60) * ϕ₆    #
  )

  # FiniteDifferenceMethod(-4:1, 1)
  ∂ϕ₅ = (
    (1 / 20) * ϕ₁ -  #
    (1 / 3) * ϕ₂ +   #
    ϕ₃ -             #
    2ϕ₄ +            #
    (13 / 12) * ϕ₅ + #
    (1 / 5) * ϕ₆     #
  )

  # FiniteDifferenceMethod(-3:2, 1)
  ∂ϕ₄ = (
    (-1 / 30) * ϕ₁ + #
    (1 / 4) * ϕ₂ -   #
    ϕ₃ +             #
    (1 / 3) * ϕ₄ +   #
    (1 / 2) * ϕ₅ -   #
    (1 / 20) * ϕ₆    #
  )

  return (∂ϕ₄, ∂ϕ₅, ∂ϕ₆)
end

"""
Second derivative of ϕ, using the first derivatives ∂ϕ and cell-center quantities ϕ
"""
function ∂²(∂ϕᵢ₊₁, ∂ϕᵢ₋₁, ϕᵢ₋₁, ϕᵢ, ϕᵢ₊₁)
  ∂²ϕᵢ = 2((ϕᵢ₊₁ + ϕᵢ₋₁) - 2ϕᵢ) - (∂ϕᵢ₊₁ - ∂ϕᵢ₋₁) / 2
  return ∂²ϕᵢ
end

@inline ϕᴸᵢ₊½(ϕᵢ, ∂ϕᵢ, ∂²ϕᵢ) = ϕᵢ + ((1 / 2) * ∂ϕᵢ + (1 / 12) * ∂²ϕᵢ)
@inline ϕᴿᵢ₋½(ϕᵢ, ∂ϕᵢ, ∂²ϕᵢ) = ϕᵢ - ((1 / 2) * ∂ϕᵢ + (1 / 12) * ∂²ϕᵢ)
@inline ϕᴿᵢ₊½(ϕᵢ₊₁, ∂ϕᵢ₊₁, ∂²ϕᵢ₊₁) = ϕᵢ₊₁ - ((1 / 2) * ∂ϕᵢ₊₁ + (1 / 12) * ∂²ϕᵢ₊₁)

function meg6_ϕᵢ₊½(ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅)
  (607 / 960) * (ϕᵢ + ϕᵢ₊₁) - (461 / 2880) * (ϕᵢ₊₂ + ϕᵢ₋₁) + (17 / 576) * (ϕᵢ₊₃ + ϕᵢ₋₂) -
  (1 / 720) * (ϕᵢ₊₄ + ϕᵢ₋₃) - (1 / 2880) * (ϕᵢ₊₅ + ϕᵢ₋₄)
end

function meg6_ϕᴸᴿᵢ₊½(ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅)
  ϕᴸᵢ₊½ = (
    +(35 / 48) * ϕᵢ                #
    + (257 / 480) * ϕᵢ₊₁           #
    - (19 / 180) * ϕᵢ₊₂            #
    + (7 / 480) * ϕᵢ₊₃             #
    - (1 / 1440) * (ϕᵢ₊₄ + ϕᵢ₋₄)   #
    - (103 / 480) * ϕᵢ₋₁           #
    + (2 / 45) * ϕᵢ₋₂              #
    - (1 / 480) * ϕᵢ₋₃             #
  )

  ϕᴿᵢ₊½ = (
    +(257 / 480) * ϕᵢ              #
    + (35 / 48) * ϕᵢ₊₁             #
    - (103 / 480) * ϕᵢ₊₂           #
    + (2 / 45) * ϕᵢ₊₃              #
    - (1 / 480) * ϕᵢ₊₄             #
    - (1 / 1440) * (ϕᵢ₊₅ + ϕᵢ₋₃)   #
    - (19 / 180) * ϕᵢ₋₁            #
    + (7 / 480) * ϕᵢ₋₂             #
  )

  return (ϕᴸᵢ₊½, ϕᴿᵢ₊½)
end

"""
    central_deriv_6th_order(ϕᵢ₋₅, ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅)

Find the first derivative at the cell-center, e.g., ∂ϕ/∂ξ|ᵢ based on a set of neighbor ϕ using a 6th order scheme. 
This can be used to compute grid metrics or other derivatives on the inner region of a grid for instance. Edges of
the grid need to use edge derivative functions
"""
function inner_central_deriv_6th_order(
  ϕᵢ₋₅, ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅
)
  ∂ϕᵢ = (
    (1141 / 1440) * (ϕᵢ₊₁ - ϕᵢ₋₁) +
    (91 / 480) * (ϕᵢ₋₂ - ϕᵢ₊₂) +
    (89 / 2880) * (ϕᵢ₊₃ - ϕᵢ₋₃) +
    (1 / 960) * (ϕᵢ₋₄ - ϕᵢ₊₄) +
    (1 / 2880) * (ϕᵢ₋₅ - ϕᵢ₊₅)
  )

  return ∂ϕᵢ
end

"""
    interface_deriv_6th_order_alpha_damping(
  ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅
)

Find the first derivative at the interface, e.g., ∂ϕ/∂ξ|ᵢ₊½ based on a set neighbor ϕ using a 6th order scheme. This 
uses alpha damping, which is designed to damp out high-frequency modes.
"""
function inner_interface_deriv_6th_order_alpha_damping(
  ϕᵢ₋₄, ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄, ϕᵢ₊₅
)
  ∂ϕᵢ₊½ = (
    (61 / 80) * (ϕᵢ₊₁ - ϕᵢ) +
    (59 / 720) * (ϕᵢ₊₂ - ϕᵢ₋₁) +
    (1 / 144) * (ϕᵢ₋₂ - ϕᵢ₊₃) +
    (1 / 180) * (ϕᵢ₊₄ - ϕᵢ₋₃) +
    (1 / 720) * (ϕᵢ₋₄ - ϕᵢ₊₅)
  )

  return ∂ϕᵢ₊½
end

"""
    interface_deriv_6th_order(ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄)

"""
function inner_interface_deriv_6th_order(ϕᵢ₋₃, ϕᵢ₋₂, ϕᵢ₋₁, ϕᵢ, ϕᵢ₊₁, ϕᵢ₊₂, ϕᵢ₊₃, ϕᵢ₊₄)
  ∂ϕᵢ₊½ = (
    (3 / 8) * (ϕᵢ₊₁ - ϕᵢ) +
    (3 / 10) * (ϕᵢ₊₂ - ϕᵢ₋₁) +
    (1 / 15) * (ϕᵢ₋₂ - ϕᵢ₊₃) +
    (1 / 120) * (ϕᵢ₊₄ - ϕᵢ₋₃)
  )

  return ∂ϕᵢ₊½
end

# @kernel inbounds = true function ∂²_kernel!(
#   ∂²W::AbstractArray{T,N}, # 2nd derivatives
#   ∂W, # 1st derivatives
#   W, # cell-centered average value
#   axis, # axis we're doing the derivative along
#   I0, # initial offset due to halo regions
# ) where {T,N}
#   I = @index(Global, Cartesian)
#   I += I0

#   ᵢ₊₁ = shift(I, axis, +1)
#   ᵢ₋₁ = shift(I, axis, -1)

#   dW_1 = ∂W[ᵢ₊₁] - ∂W[ᵢ₋₁]
#   dW_1 = dW_1 * !isapprox(W[ᵢ₊₁], W[ᵢ₋₁])

#   ∂²W[I] = 2((W[ᵢ₊₁] + W[ᵢ₋₁]) - 2W[I]) - dW_1 / 2
# end

# @kernel inbounds = true function ∂_kernel!(
#   ∂W::AbstractArray{T,N},
#   W,
#   axis, # axis we're doing the derivative along
#   I0, # initial offset due to halo regions
# ) where {T,N}
#   f1 = T(3 / 4) # convert to proper datatype, e.g. Float32 is faster on GPUs
#   f2 = T(3 / 20)
#   f3 = T(1 / 60)

#   I = @index(Global, Cartesian)
#   I += I0

#   ᵢ₊₁ = shift(I, axis, +1)
#   ᵢ₊₂ = shift(I, axis, +2)
#   ᵢ₊₃ = shift(I, axis, +3)
#   ᵢ₋₁ = shift(I, axis, -1)
#   ᵢ₋₂ = shift(I, axis, -2)
#   ᵢ₋₃ = shift(I, axis, -3)

#   dW_1 = W[ᵢ₊₁] - W[ᵢ₋₁]
#   dW_2 = W[ᵢ₊₂] - W[ᵢ₋₂]
#   dW_3 = W[ᵢ₊₃] - W[ᵢ₋₃]

#   dW_1 = dW_1 * !isapprox(W[ᵢ₊₁], W[ᵢ₋₁])
#   dW_2 = dW_2 * !isapprox(W[ᵢ₊₂], W[ᵢ₋₂])
#   dW_3 = dW_3 * !isapprox(W[ᵢ₊₃], W[ᵢ₋₃])

#   ∂W[I] = (f1 * dW_1 - f2 * dW_2 + f3 * dW_3)
# end
end