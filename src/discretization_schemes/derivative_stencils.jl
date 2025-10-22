@inline function kahan_sum(values)
  s = zero(eltype(values))
  c = zero(eltype(values))
  for t in values
    y = t - c
    tmp = s + y
    c = (tmp - s) - y
    s = tmp
  end
  return s + c
end

@inline function tol_diff(a::T, b::T) where {T}
  a_m_b = a - b
  return a_m_b * !isapprox(a, b) #; rtol=1e-7)
end

# ------------------------
# Central derivatives
# ------------------------
# 2nd order
function central_derivative(ϕᵢ₋₁::T, ϕᵢ₊₁::T) where {T}
  ∂ϕ = tol_diff(ϕᵢ₊₁, ϕᵢ₋₁) / 2
  return ∂ϕ
end

# 4th order
function central_derivative(ϕᵢ₋₂::T, ϕᵢ₋₁::T, ϕᵢ₊₁::T, ϕᵢ₊₂::T) where {T}
  ∂ϕ = kahan_sum(((2 / 3) * tol_diff(ϕᵢ₊₁, ϕᵢ₋₁), (-1 / 12) * tol_diff(ϕᵢ₊₂, ϕᵢ₋₂)))
  return ∂ϕ
end

# 6th order
function central_derivative(ϕᵢ₋₃::T, ϕᵢ₋₂::T, ϕᵢ₋₁::T, ϕᵢ₊₁::T, ϕᵢ₊₂::T, ϕᵢ₊₃::T) where {T}
  ∂ϕ = kahan_sum((
    (3 / 4) * tol_diff(ϕᵢ₊₁, ϕᵢ₋₁),   #
    (-3 / 20) * tol_diff(ϕᵢ₊₂, ϕᵢ₋₂), #
    (1 / 60) * tol_diff(ϕᵢ₊₃, ϕᵢ₋₃),  #
  ))
  return ∂ϕ
end

# 8th order
function central_derivative(
  ϕᵢ₋₄::T, ϕᵢ₋₃::T, ϕᵢ₋₂::T, ϕᵢ₋₁::T, ϕᵢ₊₁::T, ϕᵢ₊₂::T, ϕᵢ₊₃::T, ϕᵢ₊₄::T
) where {T}
  ∂ϕ = kahan_sum((
    (4 / 5) * tol_diff(ϕᵢ₊₁, ϕᵢ₋₁),     #
    (-1 / 5) * tol_diff(ϕᵢ₊₂, ϕᵢ₋₂),    #
    (4 / 105) * tol_diff(ϕᵢ₊₃, ϕᵢ₋₃),   #
    (-1 / 280) * tol_diff(ϕᵢ₊₄, ϕᵢ₋₄),  #
  ))
  return ∂ϕ
end

# ------------------------
# One-sided derivatives
# ------------------------
# 2nd order
function forward_derivative(ϕᵢ::T, ϕᵢ₊₁::T; ϵ::T=eps(T)) where {T}
  ∂ϕ = ϕᵢ₊₁ - ϕᵢ
  return ∂ϕ * (abs(∂ϕ) >= ϵ)
end

# 4th order
function forward_derivative(ϕᵢ::T, ϕᵢ₊₁::T, ϕᵢ₊₂::T; ϵ::T=eps(T)) where {T}
  ∂ϕ = kahan_sum((-3ϕᵢ, 4ϕᵢ₊₁, -ϕᵢ₊₂)) / 2
  return ∂ϕ * (abs(∂ϕ) >= ϵ)
end

# 6th order
function forward_derivative(ϕᵢ::T, ϕᵢ₊₁::T, ϕᵢ₊₂::T, ϕᵢ₊₃::T; ϵ::T=eps(T)) where {T}
  ∂ϕ = kahan_sum((-11ϕᵢ, 18ϕᵢ₊₁, -9ϕᵢ₊₂, 2ϕᵢ₊₃)) / 6
  return ∂ϕ * (abs(∂ϕ) >= ϵ)
end

# 8th order
function forward_derivative(
  ϕᵢ::T, ϕᵢ₊₁::T, ϕᵢ₊₂::T, ϕᵢ₊₃::T, ϕᵢ₊₄::T; ϵ::T=eps(T)
) where {T}
  ∂ϕ = kahan_sum((-25ϕᵢ, 48ϕᵢ₊₁, -36ϕᵢ₊₂, 16ϕᵢ₊₃, -3ϕᵢ₊₄)) / 12
  return ∂ϕ * (abs(∂ϕ) >= ϵ)
end

function backward_derivative(args...)
  ∂ϕ = -forward_derivative(reverse(args)...)
  return ∂ϕ
end

# ------------------------
# Mixed derivatives
# ------------------------
function mixed_derivatives_loedge(
  ϕ₁::T, ϕ₂::T, ϕ₃::T, ϕ₄::T, ϕ₅::T, ϕ₆::T; ϵ::T=50eps(T)
) where {T}
  # FiniteDifferenceMethod(0:5, 1)
  ∂ϕ₁ = kahan_sum((
    (-137 / 60) * ϕ₁,  #
    5ϕ₂,               #
    -5ϕ₃,              #
    (10 / 3) * ϕ₄,     #
    -(5 / 4) * ϕ₅,     #
    (1 / 5) * ϕ₆,      #
  ))

  # FiniteDifferenceMethod(-1:4, 1)
  ∂ϕ₂ = kahan_sum((
    (-1 / 5) * ϕ₁,    #
    -(13 / 12) * ϕ₂,  #
    2ϕ₃,              #
    -ϕ₄,              #
    (1 / 3) * ϕ₅,     #
    -(1 / 20) * ϕ₆,   #
  ))

  # FiniteDifferenceMethod(-2:3, 1)
  ∂ϕ₃ = kahan_sum((
    (1 / 20) * ϕ₁,  # 
    -(1 / 2) * ϕ₂,  #
    -(1 / 3) * ϕ₃,  #
    ϕ₄,             #
    -(1 / 4) * ϕ₅,  #
    (1 / 30) * ϕ₆,  #
  ))

  return (
    ∂ϕ₁ * (abs(∂ϕ₁) >= ϵ),#
    ∂ϕ₂ * (abs(∂ϕ₂) >= ϵ),#
    ∂ϕ₃ * (abs(∂ϕ₃) >= ϵ), #
  )
end

function mixed_derivatives_hiedge(
  ϕ₁::T, ϕ₂::T, ϕ₃::T, ϕ₄::T, ϕ₅::T, ϕ₆::T; ϵ::T=50eps(T)
) where {T}
  # FiniteDifferenceMethod(-5:0, 1)
  ∂ϕ₆ = kahan_sum((
    (-1 / 5) * ϕ₁,   #
    (5 / 4) * ϕ₂,    #
    -(10 / 3) * ϕ₃,  #
    5ϕ₄,             #
    -5ϕ₅,            #
    (137//60) * ϕ₆,  #
  ))

  # FiniteDifferenceMethod(-4:1, 1)
  ∂ϕ₅ = kahan_sum((
    (1 / 20) * ϕ₁,   #
    -(1 / 3) * ϕ₂,   #
    ϕ₃,              #
    -2ϕ₄,            #
    (13 / 12) * ϕ₅,  #
    (1 / 5) * ϕ₆,    #
  ))

  # FiniteDifferenceMethod(-3:2, 1)
  ∂ϕ₄ = kahan_sum((
    (-1 / 30) * ϕ₁,  #
    (1 / 4) * ϕ₂,    #
    -ϕ₃,             #
    (1 / 3) * ϕ₄,    #
    (1 / 2) * ϕ₅,    #
    -(1 / 20) * ϕ₆,  #
  ))

  return (
    ∂ϕ₄ * (abs(∂ϕ₄) >= ϵ),#
    ∂ϕ₅ * (abs(∂ϕ₅) >= ϵ),#
    ∂ϕ₆ * (abs(∂ϕ₆) >= ϵ), #
  )
end

# ------------------------
# Second derivative
# ------------------------
function second_derivative(ϕᵢ₋₁::T, ϕᵢ::T, ϕᵢ₊₁::T, ∂ϕᵢ₋₁, ∂ϕᵢ₊₁; ϵ::T=50eps(T)) where {T}
  # ∂²ϕ = 2((ϕᵢ₊₁ + ϕᵢ₋₁) - 2ϕᵢ) - (∂ϕᵢ₊₁ - ∂ϕᵢ₋₁) / 2

  d_∂ϕ = -tol_diff(∂ϕᵢ₊₁, ∂ϕᵢ₋₁) / 2
  # ∂²ϕ = 2((ϕᵢ₊₁ + ϕᵢ₋₁) - 2ϕᵢ) - dϕ / 2

  ∂²ϕ = kahan_sum((2(ϕᵢ₊₁ + ϕᵢ₋₁), -4ϕᵢ, d_∂ϕ))
  # ∂²ϕ = kahan_sum((2ϕᵢ₊₁, 2ϕᵢ₋₁, -4ϕᵢ, -∂ϕᵢ₊₁ / 2, ∂ϕᵢ₋₁ / 2))
  # if 0 < abs(∂²ϕ) < eps()
  #   @show ϕᵢ₋₁ ϕᵢ ϕᵢ₊₁ ∂ϕᵢ₋₁ ∂ϕᵢ₊₁ d_∂ϕ ∂²ϕ
  #   error("what the bork")
  # end
  return ∂²ϕ #* (abs(∂²ϕ) >= ϵ)
end
