#-------------------------------------------------------------
# 3D cell-center derivatives
#-------------------------------------------------------------
∂ϕ_∂ξ_3d(ϕ, backend) = ∂ϕ_∂ξ_3d(ϕ, backend, EdgeInterpolationOrder3())

function ∂ϕ_∂ξ_3d(ϕ, backend, ::EdgeInterpolationOrder1)
  ϕᵢ₊½(t, i, j, k, p) = _edge_reconstruct(
    ϕ(t, i, j, k, p), ϕ(t, i + 1, j, k, p), EdgeInterpolationOrder1()
  )
  ∂ϕ_∂ξ(t, i, j, k, p) = ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  return ∂ϕ_∂ξ
end

function ∂ϕ_∂ξ_3d(ϕ, backend, ::EdgeInterpolationOrder2)
  ξ_derivs(t, i, j, k, p) = value_and_derivative(ξ -> ϕ(t, ξ, j, k, p), backend, i)

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    return _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, EdgeInterpolationOrder2())
  end

  ∂ϕ_∂ξ(t, i, j, k, p) = ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  return ∂ϕ_∂ξ
end

function ∂ϕ_∂ξ_3d(ϕ, backend, ::EdgeInterpolationOrder3)
  ξ_derivs(t, i, j, k, p) = value_derivative_and_second_derivative(
    ξ -> ϕ(t, ξ, j, k, p), backend, i
  )

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    return _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, EdgeInterpolationOrder3()
    )
  end

  ∂ϕ_∂ξ(t, i, j, k, p) = ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  return ∂ϕ_∂ξ
end

∂ϕ_∂η_3d(ϕ, backend) = ∂ϕ_∂η_3d(ϕ, backend, EdgeInterpolationOrder3())

function ∂ϕ_∂η_3d(ϕ, backend, ::EdgeInterpolationOrder1)
  ϕⱼ₊½(t, i, j, k, p) = _edge_reconstruct(
    ϕ(t, i, j, k, p), ϕ(t, i, j + 1, k, p), EdgeInterpolationOrder1()
  )
  ∂ϕ_∂η(t, i, j, k, p) = ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  return ∂ϕ_∂η
end

function ∂ϕ_∂η_3d(ϕ, backend, ::EdgeInterpolationOrder2)
  η_derivs(t, i, j, k, p) = value_and_derivative(η -> ϕ(t, i, η, k, p), backend, j)

  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    return _edge_reconstruct(ϕⱼ, ∂ϕ_∂ηⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, EdgeInterpolationOrder2())
  end

  ∂ϕ_∂η(t, i, j, k, p) = ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  return ∂ϕ_∂η
end

function ∂ϕ_∂η_3d(ϕ, backend, ::EdgeInterpolationOrder3)
  η_derivs(t, i, j, k, p) = value_derivative_and_second_derivative(
    η -> ϕ(t, i, η, k, p), backend, j
  )

  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    return _edge_reconstruct(
      ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁, EdgeInterpolationOrder3()
    )
  end

  ∂ϕ_∂η(t, i, j, k, p) = ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  return ∂ϕ_∂η
end

∂ϕ_∂ζ_3d(ϕ, backend) = ∂ϕ_∂ζ_3d(ϕ, backend, EdgeInterpolationOrder3())

function ∂ϕ_∂ζ_3d(ϕ, backend, ::EdgeInterpolationOrder1)
  ϕₖ₊½(t, i, j, k, p) = _edge_reconstruct(
    ϕ(t, i, j, k, p), ϕ(t, i, j, k + 1, p), EdgeInterpolationOrder1()
  )
  ∂ϕ_∂ζ(t, i, j, k, p) = ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)
  return ∂ϕ_∂ζ
end

function ∂ϕ_∂ζ_3d(ϕ, backend, ::EdgeInterpolationOrder2)
  ζ_derivs(t, i, j, k, p) = value_and_derivative(ζ -> ϕ(t, i, j, ζ, p), backend, k)

  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    return _edge_reconstruct(ϕₖ, ∂ϕ_∂ζₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, EdgeInterpolationOrder2())
  end

  ∂ϕ_∂ζ(t, i, j, k, p) = ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)
  return ∂ϕ_∂ζ
end

function ∂ϕ_∂ζ_3d(ϕ, backend, ::EdgeInterpolationOrder3)
  ζ_derivs(t, i, j, k, p) = value_derivative_and_second_derivative(
    ζ -> ϕ(t, i, j, ζ, p), backend, k
  )

  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    return _edge_reconstruct(
      ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁, EdgeInterpolationOrder3()
    )
  end

  ∂ϕ_∂ζ(t, i, j, k, p) = ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)
  return ∂ϕ_∂ζ
end
