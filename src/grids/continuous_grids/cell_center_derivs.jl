
#-------------------------------------------------------------
# 1D cell-center derivatives
#-------------------------------------------------------------
function ∂ϕ_∂ξ_1d(ϕ, backend)
  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = value_derivative_and_second_derivative(
      ξ -> ϕ(t, ξ, j, k, p), backend, i
    )
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = value_derivative_and_second_derivative(
      ξ -> ϕ(t, ξ, j, k, p), backend, i + 1
    )

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ∂ϕ_∂ξ(t, i, j, k, p)
    ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  end

  return ∂ϕ_∂ξ
end

#-------------------------------------------------------------
# 2D cell-center derivatives
#-------------------------------------------------------------
function ∂ϕ_∂ξ_2d(ϕ, backend)
  function ϕᵢ₊½(t::T, i::I, j::J, p) where {T,I,J}
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = value_derivative_and_second_derivative(
      ξ -> ϕ(t, ξ, j, p), backend, i
    )
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = value_derivative_and_second_derivative(
      ξ -> ϕ(t, ξ, j, p), backend, i + 1
    )

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ∂ϕ_∂ξ(t::T, i::I, j::J, p) where {T,I,J}
    ϕᵢ₊½(t, i, j, p) - ϕᵢ₊½(t, i - 1, j, p)
  end

  return ∂ϕ_∂ξ
end

function ∂ϕ_∂η_2d(ϕ, backend)
  function ϕⱼ₊½(t::T, i::I, j::J, p) where {T,I,J}
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = value_derivative_and_second_derivative(
      η -> ϕ(t, i, η, p), backend, j
    )
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = value_derivative_and_second_derivative(
      η -> ϕ(t, i, η, p), backend, j + 1
    )

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  function ∂ϕ_∂η(t, i, j, p)
    ϕⱼ₊½(t, i, j, p) - ϕⱼ₊½(t, i, j - 1, p)
  end

  return ∂ϕ_∂η
end

#-------------------------------------------------------------
# 3D cell-center derivatives
#-------------------------------------------------------------
function ∂ϕ_∂ξ_3d(ϕ, backend)
  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = value_derivative_and_second_derivative(
      ξ -> ϕ(t, ξ, j, k, p), backend, i
    )
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = value_derivative_and_second_derivative(
      ξ -> ϕ(t, ξ, j, k, p), backend, i + 1
    )

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ∂ϕ_∂ξ(t, i, j, k, p)
    ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  end

  return ∂ϕ_∂ξ
end

function ∂ϕ_∂η_3d(ϕ, backend)
  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = value_derivative_and_second_derivative(
      η -> ϕ(t, i, η, k, p), backend, j
    )
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = value_derivative_and_second_derivative(
      η -> ϕ(t, i, η, k, p), backend, j + 1
    )

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  function ∂ϕ_∂η(t, i, j, k, p)
    ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  end

  return ∂ϕ_∂η
end

function ∂ϕ_∂ζ_3d(ϕ, backend)
  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ξₖ, ∂²ϕ_∂ξ²ₖ = value_derivative_and_second_derivative(
      ζ -> ϕ(t, i, j, ζ, p), backend, k
    )
    ϕₖ₊₁, ∂ϕ_∂ξₖ₊₁, ∂²ϕ_∂ξ²ₖ₊₁ = value_derivative_and_second_derivative(
      ζ -> ϕ(t, i, j, ζ, p), backend, k + 1
    )

    ϕᴸₖ₊½ = ϕₖ + (1 / 2) * ∂ϕ_∂ξₖ + (1 / 12) * ∂²ϕ_∂ξ²ₖ
    ϕᴿₖ₊½ = ϕₖ₊₁ - (1 / 2) * ∂ϕ_∂ξₖ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ₖ₊₁

    return (ϕᴸₖ₊½ + ϕᴿₖ₊½) / 2
  end

  function ∂ϕ_∂ζ(t, i, j, k, p)
    ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)
  end

  return ∂ϕ_∂ζ
end