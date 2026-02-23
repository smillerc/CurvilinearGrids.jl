struct MetricCache{FM,IM,EM}
  forward::FM
  inverse::IM
  edge::EM
end

abstract type EdgeInterpolationSchemeTrait end

struct EdgeInterpolationOrder1 <: EdgeInterpolationSchemeTrait end
struct EdgeInterpolationOrder2 <: EdgeInterpolationSchemeTrait end
struct EdgeInterpolationOrder3 <: EdgeInterpolationSchemeTrait end

@inline _edge_reconstruct(ϕᵢ, ϕᵢ₊₁, ::EdgeInterpolationOrder1) = 0.5 * (ϕᵢ + ϕᵢ₊₁)

@inline function _edge_reconstruct(
  ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ::EdgeInterpolationOrder2, Δξ::Real=1
)
  h = 0.5 * Δξ
  ϕᴸᵢ₊½ = ϕᵢ + h * ∂ϕ_∂ξᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - h * ∂ϕ_∂ξᵢ₊₁
  return 0.5 * (ϕᴸᵢ₊½ + ϕᴿᵢ₊½)
end

@inline function _edge_reconstruct(
  ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, ::EdgeInterpolationOrder3, Δξ::Real=1
)
  h = 0.5 * Δξ
  ϕᴸᵢ₊½ = ϕᵢ + h * ∂ϕ_∂ξᵢ + (Δξ^2 / 12) * ∂²ϕ_∂ξ²ᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - h * ∂ϕ_∂ξᵢ₊₁ + (Δξ^2 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁
  return 0.5 * (ϕᴸᵢ₊½ + ϕᴿᵢ₊½)
end

"""
3D metric cache
"""
function MetricCache(
  x::Function,
  y::Function,
  z::Function,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  xξ(t, i, j, k, p) = derivative(ξ -> x(t, ξ, j, k, p), backend, i)
  xη(t, i, j, k, p) = derivative(η -> x(t, i, η, k, p), backend, j)
  xζ(t, i, j, k, p) = derivative(ζ -> x(t, i, j, ζ, p), backend, k)
  yξ(t, i, j, k, p) = derivative(ξ -> y(t, ξ, j, k, p), backend, i)
  yη(t, i, j, k, p) = derivative(η -> y(t, i, η, k, p), backend, j)
  yζ(t, i, j, k, p) = derivative(ζ -> y(t, i, j, ζ, p), backend, k)
  zξ(t, i, j, k, p) = derivative(ξ -> z(t, ξ, j, k, p), backend, i)
  zη(t, i, j, k, p) = derivative(η -> z(t, i, η, k, p), backend, j)
  zζ(t, i, j, k, p) = derivative(ζ -> z(t, i, j, ζ, p), backend, k)

  function jacobian_matrix(t, i, j, k, p)

    # compute the jacobian matrix w/o any extra logic
    function compute_jacobian_matrix(t, i, j, k, p)
      jac = DifferentiationInterface.jacobian(
        u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
        backend,
        @SVector [i, j, k]
      )

      return jac
    end

    jac_matrix = compute_jacobian_matrix(t, i, j, k, p)
    return jac_matrix
    # J = det(jac_matrix)

  end

  forward_metrics = (;
    jacobian=jacobian_matrix,
    J=jacobian,
    xξ=xξ,
    xη=xη,
    xζ=xζ,
    yξ=yξ,
    yη=yη,
    yζ=yζ,
    zξ=zξ,
    zη=zη,
    zζ=zζ,
  )

  inverse_metrics, edge_metrics = get_inverse_metric_terms(
    x, y, z, backend; edge_interpolation_scheme=edge_interpolation_scheme
  )

  # edge_metrics = get_edge_functions_3d(forward_metrics, inverse_metrics, backend)

  return MetricCache(forward_metrics, inverse_metrics, edge_metrics)
end

"""
2D metric cache
"""
function MetricCache(
  x::Function,
  y::Function,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  xξ(t, i, j, p) = derivative(ξ -> x(t, ξ, j, p), backend, i)
  xη(t, i, j, p) = derivative(η -> x(t, i, η, p), backend, j)
  xτ(t, i, j, p) = derivative(τ -> x(τ, i, j, p), backend, t)

  yξ(t, i, j, p) = derivative(ξ -> y(t, ξ, j, p), backend, i)
  yη(t, i, j, p) = derivative(η -> y(t, i, η, p), backend, j)
  yτ(t, i, j, p) = derivative(τ -> y(τ, i, j, p), backend, t)

  function jacobian_matrix(t, i, j, p)

    # compute the jacobian matrix w/o any extra logic
    # function compute_jacobian_matrix(t, i, j, p)
    jac = DifferentiationInterface.jacobian(
      u -> SVector(x(t, u..., p), y(t, u..., p)), backend, @SVector [i, j]
    )

    return jac
  end

  jacobian(t, i, j, p) = det(jacobian_matrix(t, i, j, p))
  function jinv(t, i, j, p)
    J_inv = inv(jacobian_matrix(t, i, j, p))
    return J_inv
  end

  function normalized_jinv(t, i, j, p)
    jac = jacobian_matrix(t, i, j, p)
    J = det(jac)
    J_inv = inv(jac)

    return J_inv .* J
  end

  forward_metrics = (; jacobian=jacobian_matrix, J=jacobian, xξ=xξ, xη=xη, yξ=yξ, yη=yη)

  inverse_metrics = (;
    # ξ̂x=ξ̂x, ξ̂y=ξ̂y, η̂x=η̂x, η̂y=η̂y,
    Jinv=jinv,
    Jinv_norm=normalized_jinv,
  )

  edge_metrics = get_edge_functions_2d(
    forward_metrics,
    inverse_metrics,
    backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )

  return MetricCache(forward_metrics, inverse_metrics, edge_metrics)
end

"""
1D metric cache
"""
function MetricCache(
  x::Function,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  xξ(t, i, p) = derivative(ξ -> x(t, ξ, p), backend, i)

  jacobian_matrix(t, i, p) = @SMatrix [xξ(t, i, p)]

  jacobian(t, i, p) = det(jacobian_matrix(t, i, p))
  jinv(t, i, p) = inv(jacobian_matrix(t, i, p))

  function normalized_jinv(t, i, p)
    jac = jacobian_matrix(t, i, p)
    J = det(jac)
    J_inv = inv(jac)

    return J_inv .* J
  end

  forward_metrics = (; jacobian=jacobian_matrix, J=jacobian, xξ=xξ)

  inverse_metrics = (;
    # ξ̂x=ξ̂x, ξ̂y=ξ̂y, η̂x=η̂x, η̂y=η̂y,
    Jinv=jinv,
    Jinv_norm=normalized_jinv,
  )

  edge_metrics = get_edge_functions_1d(
    forward_metrics,
    inverse_metrics,
    backend;
    edge_interpolation_scheme=edge_interpolation_scheme,
  )

  return MetricCache(forward_metrics, inverse_metrics, edge_metrics)
end

edge_functions_3d(ϕ, backend) = edge_functions_3d(ϕ, backend, EdgeInterpolationOrder3())

function edge_functions_3d(ϕ, backend, ::EdgeInterpolationOrder1)
  function ϕᵢ₊½(t, i, j, k, p)
    _edge_reconstruct(ϕ(t, i, j, k, p), ϕ(t, i + 1, j, k, p), EdgeInterpolationOrder1())
  end
  function ϕⱼ₊½(t, i, j, k, p)
    _edge_reconstruct(ϕ(t, i, j, k, p), ϕ(t, i, j + 1, k, p), EdgeInterpolationOrder1())
  end
  function ϕₖ₊½(t, i, j, k, p)
    _edge_reconstruct(ϕ(t, i, j, k, p), ϕ(t, i, j, k + 1, p), EdgeInterpolationOrder1())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

function edge_functions_3d(ϕ, backend, ::EdgeInterpolationOrder2)
  ξ_derivs(t, i, j, k, p) = value_and_derivative(ξ -> ϕ(t, ξ, j, k, p), backend, i)
  η_derivs(t, i, j, k, p) = value_and_derivative(η -> ϕ(t, i, η, k, p), backend, j)
  ζ_derivs(t, i, j, k, p) = value_and_derivative(ζ -> ϕ(t, i, j, ζ, p), backend, k)

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, EdgeInterpolationOrder2())
  end
  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    _edge_reconstruct(ϕⱼ, ∂ϕ_∂ηⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, EdgeInterpolationOrder2())
  end
  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    _edge_reconstruct(ϕₖ, ∂ϕ_∂ζₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, EdgeInterpolationOrder2())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

function edge_functions_3d(ϕ, backend, ::EdgeInterpolationOrder3)
  function ξ_derivs(t, i, j, k, p)
    value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, j, k, p), backend, i)
  end
  function η_derivs(t, i, j, k, p)
    value_derivative_and_second_derivative(η -> ϕ(t, i, η, k, p), backend, j)
  end
  function ζ_derivs(t, i, j, k, p)
    value_derivative_and_second_derivative(ζ -> ϕ(t, i, j, ζ, p), backend, k)
  end

  function ϕᵢ₊½(t, i, j, k, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, j, k, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, j, k, p)
    _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, EdgeInterpolationOrder3()
    )
  end
  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁ = η_derivs(t, i, j + 1, k, p)
    _edge_reconstruct(
      ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁, EdgeInterpolationOrder3()
    )
  end
  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)
    _edge_reconstruct(
      ϕₖ, ∂ϕ_∂ζₖ, ∂²ϕ_∂ζ²ₖ, ϕₖ₊₁, ∂ϕ_∂ζₖ₊₁, ∂²ϕ_∂ζ²ₖ₊₁, EdgeInterpolationOrder3()
    )
  end
  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

edge_functions_2d(ϕ, backend) = edge_functions_2d(ϕ, backend, EdgeInterpolationOrder3())

function edge_functions_2d(ϕ, backend, ::EdgeInterpolationOrder1)
  function ϕᵢ₊½(t, i, j, p)
    _edge_reconstruct(ϕ(t, i, j, p), ϕ(t, i + 1, j, p), EdgeInterpolationOrder1())
  end
  function ϕⱼ₊½(t, i, j, p)
    _edge_reconstruct(ϕ(t, i, j, p), ϕ(t, i, j + 1, p), EdgeInterpolationOrder1())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½)
end

function edge_functions_2d(ϕ, backend, ::EdgeInterpolationOrder2)
  ξ_derivs(t, i, j, p) = value_and_derivative(ξ -> ϕ(t, ξ, j, p), backend, i)
  η_derivs(t, i, j, p) = value_and_derivative(η -> ϕ(t, i, η, p), backend, j)

  function ϕᵢ₊½(t, i, j, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, j, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, j, p)
    _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, EdgeInterpolationOrder2())
  end
  function ϕⱼ₊½(t, i, j, p)
    ϕⱼ, ∂ϕ_∂ηⱼ = η_derivs(t, i, j, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁ = η_derivs(t, i, j + 1, p)
    _edge_reconstruct(ϕⱼ, ∂ϕ_∂ηⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, EdgeInterpolationOrder2())
  end
  return (; ϕᵢ₊½, ϕⱼ₊½)
end

function edge_functions_2d(ϕ, backend, ::EdgeInterpolationOrder3)
  function ξ_derivs(t, i, j, p)
    value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, j, p), backend, i)
  end
  function η_derivs(t, i, j, p)
    value_derivative_and_second_derivative(η -> ϕ(t, i, η, p), backend, j)
  end

  function ϕᵢ₊½(t, i, j, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, j, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, j, p)
    _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, EdgeInterpolationOrder3()
    )
  end
  function ϕⱼ₊½(t, i, j, p)
    ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ = η_derivs(t, i, j, p)
    ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁ = η_derivs(t, i, j + 1, p)
    _edge_reconstruct(
      ϕⱼ, ∂ϕ_∂ηⱼ, ∂²ϕ_∂η²ⱼ, ϕⱼ₊₁, ∂ϕ_∂ηⱼ₊₁, ∂²ϕ_∂η²ⱼ₊₁, EdgeInterpolationOrder3()
    )
  end
  return (; ϕᵢ₊½, ϕⱼ₊½)
end

edge_functions_1d(ϕ, backend) = edge_functions_1d(ϕ, backend, EdgeInterpolationOrder3())

function edge_functions_1d(ϕ, backend, ::EdgeInterpolationOrder1)
  ϕᵢ₊½(t, i, p) = _edge_reconstruct(ϕ(t, i, p), ϕ(t, i + 1, p), EdgeInterpolationOrder1())
  return (; ϕᵢ₊½)
end

function edge_functions_1d(ϕ, backend, ::EdgeInterpolationOrder2)
  ξ_derivs(t, i, p) = value_and_derivative(ξ -> ϕ(t, ξ, p), backend, i)
  function ϕᵢ₊½(t, i, p)
    ϕᵢ, ∂ϕ_∂ξᵢ = ξ_derivs(t, i, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁ = ξ_derivs(t, i + 1, p)
    _edge_reconstruct(ϕᵢ, ∂ϕ_∂ξᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, EdgeInterpolationOrder2())
  end
  return (; ϕᵢ₊½)
end

function edge_functions_1d(ϕ, backend, ::EdgeInterpolationOrder3)
  ξ_derivs(t, i, p) = value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, p), backend, i)
  function ϕᵢ₊½(t, i, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, p)
    _edge_reconstruct(
      ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ, ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁, EdgeInterpolationOrder3()
    )
  end
  return (; ϕᵢ₊½)
end

function get_edge_functions_3d(
  forward_metrics,
  inverse_metrics,
  diff_backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  ξ̂xᵢ₊½, ξ̂xⱼ₊½, ξ̂xₖ₊½ = edge_functions_3d(
    inverse_metrics.ξ̂x, diff_backend, edge_interpolation_scheme
  )
  η̂xᵢ₊½, η̂xⱼ₊½, η̂xₖ₊½ = edge_functions_3d(
    inverse_metrics.η̂x, diff_backend, edge_interpolation_scheme
  )
  ζ̂xᵢ₊½, ζ̂xⱼ₊½, ζ̂xₖ₊½ = edge_functions_3d(
    inverse_metrics.ζ̂x, diff_backend, edge_interpolation_scheme
  )
  ξ̂yᵢ₊½, ξ̂yⱼ₊½, ξ̂yₖ₊½ = edge_functions_3d(
    inverse_metrics.ξ̂y, diff_backend, edge_interpolation_scheme
  )
  η̂yᵢ₊½, η̂yⱼ₊½, η̂yₖ₊½ = edge_functions_3d(
    inverse_metrics.η̂y, diff_backend, edge_interpolation_scheme
  )
  ζ̂yᵢ₊½, ζ̂yⱼ₊½, ζ̂yₖ₊½ = edge_functions_3d(
    inverse_metrics.ζ̂y, diff_backend, edge_interpolation_scheme
  )
  ξ̂zᵢ₊½, ξ̂zⱼ₊½, ξ̂zₖ₊½ = edge_functions_3d(
    inverse_metrics.ξ̂z, diff_backend, edge_interpolation_scheme
  )
  η̂zᵢ₊½, η̂zⱼ₊½, η̂zₖ₊½ = edge_functions_3d(
    inverse_metrics.η̂z, diff_backend, edge_interpolation_scheme
  )
  ζ̂zᵢ₊½, ζ̂zⱼ₊½, ζ̂zₖ₊½ = edge_functions_3d(
    inverse_metrics.ζ̂z, diff_backend, edge_interpolation_scheme
  )
  Jᵢ₊½, Jⱼ₊½, Jₖ₊½ = edge_functions_3d(
    forward_metrics.J, diff_backend, edge_interpolation_scheme
  )
  Jinv_ᵢ₊½, Jinv_ⱼ₊½, Jinv_ₖ₊½ = edge_functions_3d(
    inverse_metrics.Jinv, diff_backend, edge_interpolation_scheme
  )

  #! format: off
  edge_funcs = (;
    ξ̂xᵢ₊½, ξ̂xⱼ₊½, ξ̂xₖ₊½,         
    η̂xᵢ₊½, η̂xⱼ₊½, η̂xₖ₊½,         
    ζ̂xᵢ₊½, ζ̂xⱼ₊½, ζ̂xₖ₊½,         
    ξ̂yᵢ₊½, ξ̂yⱼ₊½, ξ̂yₖ₊½,         
    η̂yᵢ₊½, η̂yⱼ₊½, η̂yₖ₊½,         
    ζ̂yᵢ₊½, ζ̂yⱼ₊½, ζ̂yₖ₊½,         
    ξ̂zᵢ₊½, ξ̂zⱼ₊½, ξ̂zₖ₊½,         
    η̂zᵢ₊½, η̂zⱼ₊½, η̂zₖ₊½,         
    ζ̂zᵢ₊½, ζ̂zⱼ₊½, ζ̂zₖ₊½,         
    Jinv_ᵢ₊½, Jinv_ⱼ₊½, Jinv_ₖ₊½,
    Jᵢ₊½, Jⱼ₊½, Jₖ₊½,            
  )
   #! format: on

  return edge_funcs
end

function get_edge_functions_2d(
  forward_metrics,
  inverse_metrics,
  diff_backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  Jinv_ᵢ₊½, Jinv_ⱼ₊½ = edge_functions_2d(
    inverse_metrics.Jinv, diff_backend, edge_interpolation_scheme
  )
  norm_Jinv_ᵢ₊½, norm_Jinv_ⱼ₊½ = edge_functions_2d(
    inverse_metrics.Jinv_norm, diff_backend, edge_interpolation_scheme
  )

  edge_funcs = (; Jinv_ᵢ₊½, Jinv_ⱼ₊½, norm_Jinv_ᵢ₊½, norm_Jinv_ⱼ₊½)

  return edge_funcs
end

function get_edge_functions_1d(
  forward_metrics,
  inverse_metrics,
  diff_backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  Jinv_ᵢ₊½ = edge_functions_1d(
    inverse_metrics.Jinv, diff_backend, edge_interpolation_scheme
  )
  norm_Jinv_ᵢ₊½ = edge_functions_1d(
    inverse_metrics.Jinv_norm, diff_backend, edge_interpolation_scheme
  )

  edge_funcs = (; Jinv_ᵢ₊½, norm_Jinv_ᵢ₊½)

  return edge_funcs
end

function get_inverse_metric_terms(
  x,
  y,
  z,
  backend;
  edge_interpolation_scheme::EdgeInterpolationSchemeTrait=EdgeInterpolationOrder3(),
)
  edge_scheme = edge_interpolation_scheme

  #

  function get_jacobian_matrix(x, y, z, backend)
    function jacobian_matrix(t, i, j, k, p)
      DifferentiationInterface.jacobian(
        u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
        backend,
        @SVector [i, j, k]
      )
    end

    return jacobian_matrix
  end

  function get_inverse_jacobian_matrix(x, y, z, backend)
    function inverse_jacobian_matrix(t, i, j, k, p)
      return inv(
        DifferentiationInterface.jacobian(
          u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
          backend,
          @SVector [i, j, k]
        ),
      )
    end
  end

  function get_normalized_inverse_jacobian_matrix(x, y, z, backend)
    function normalized_inverse_jacobian_matrix(t, i, j, k, p)
      jac = DifferentiationInterface.jacobian(
        u -> SVector(x(t, u..., p), y(t, u..., p), z(t, u..., p)),
        backend,
        @SVector [i, j, k]
      )
      return inv(jac) .* det(jac)
    end
  end

  jacobian_matrix = get_jacobian_matrix(x, y, z, backend)
  inverse_jacobian_matrix = get_inverse_jacobian_matrix(x, y, z, backend)
  normalized_inverse_jacobian_matrix = get_normalized_inverse_jacobian_matrix(
    x, y, z, backend
  )

  xξ(t, i, j, k, p) = derivative(ξ -> x(t, ξ, j, k, p), backend, i)
  xη(t, i, j, k, p) = derivative(η -> x(t, i, η, k, p), backend, j)
  xζ(t, i, j, k, p) = derivative(ζ -> x(t, i, j, ζ, p), backend, k)
  yξ(t, i, j, k, p) = derivative(ξ -> y(t, ξ, j, k, p), backend, i)
  yη(t, i, j, k, p) = derivative(η -> y(t, i, η, k, p), backend, j)
  yζ(t, i, j, k, p) = derivative(ζ -> y(t, i, j, ζ, p), backend, k)
  zξ(t, i, j, k, p) = derivative(ξ -> z(t, ξ, j, k, p), backend, i)
  zη(t, i, j, k, p) = derivative(η -> z(t, i, η, k, p), backend, j)
  zζ(t, i, j, k, p) = derivative(ζ -> z(t, i, j, ζ, p), backend, k)

  x_ξ_y(t, i, j, k, p) = xξ(t, i, j, k, p) * y(t, i, j, k, p)
  x_η_y(t, i, j, k, p) = xη(t, i, j, k, p) * y(t, i, j, k, p)
  x_ζ_y(t, i, j, k, p) = xζ(t, i, j, k, p) * y(t, i, j, k, p)
  y_ξ_z(t, i, j, k, p) = yξ(t, i, j, k, p) * z(t, i, j, k, p)
  y_η_z(t, i, j, k, p) = yη(t, i, j, k, p) * z(t, i, j, k, p)
  y_ζ_z(t, i, j, k, p) = yζ(t, i, j, k, p) * z(t, i, j, k, p)
  z_ξ_x(t, i, j, k, p) = zξ(t, i, j, k, p) * x(t, i, j, k, p)
  z_η_x(t, i, j, k, p) = zη(t, i, j, k, p) * x(t, i, j, k, p)
  z_ζ_x(t, i, j, k, p) = zζ(t, i, j, k, p) * x(t, i, j, k, p)

  y_η_z_ζ = ∂ϕ_∂ζ_3d(y_η_z, backend, edge_scheme)
  y_ζ_z_η = ∂ϕ_∂η_3d(y_ζ_z, backend, edge_scheme)

  y_ζ_z_ξ = ∂ϕ_∂ξ_3d(y_ζ_z, backend, edge_scheme)
  y_ξ_z_ζ = ∂ϕ_∂ζ_3d(y_ξ_z, backend, edge_scheme)

  y_ξ_z_η = ∂ϕ_∂η_3d(y_ξ_z, backend, edge_scheme)
  y_η_z_ξ = ∂ϕ_∂ξ_3d(y_η_z, backend, edge_scheme)

  z_η_x_ζ = ∂ϕ_∂ζ_3d(z_η_x, backend, edge_scheme)
  z_ζ_x_η = ∂ϕ_∂η_3d(z_ζ_x, backend, edge_scheme)

  z_ζ_x_ξ = ∂ϕ_∂ξ_3d(z_ζ_x, backend, edge_scheme)
  z_ξ_x_ζ = ∂ϕ_∂ζ_3d(z_ξ_x, backend, edge_scheme)

  z_ξ_x_η = ∂ϕ_∂η_3d(z_ξ_x, backend, edge_scheme)
  z_η_x_ξ = ∂ϕ_∂ξ_3d(z_η_x, backend, edge_scheme)

  x_η_y_ζ = ∂ϕ_∂ζ_3d(x_η_y, backend, edge_scheme)
  x_ζ_y_η = ∂ϕ_∂η_3d(x_ζ_y, backend, edge_scheme)

  x_ζ_y_ξ = ∂ϕ_∂ξ_3d(x_ζ_y, backend, edge_scheme)
  x_ξ_y_ζ = ∂ϕ_∂ζ_3d(x_ξ_y, backend, edge_scheme)

  x_ξ_y_η = ∂ϕ_∂η_3d(x_ξ_y, backend, edge_scheme)
  x_η_y_ξ = ∂ϕ_∂ξ_3d(x_η_y, backend, edge_scheme)

  # Do NOT put eps() tolerance checks on these! It will create GCL-related errors
  ξ̂x(t, i, j, k, p) = y_η_z_ζ(t, i, j, k, p) − y_ζ_z_η(t, i, j, k, p)
  η̂x(t, i, j, k, p) = y_ζ_z_ξ(t, i, j, k, p) − y_ξ_z_ζ(t, i, j, k, p)
  ζ̂x(t, i, j, k, p) = y_ξ_z_η(t, i, j, k, p) − y_η_z_ξ(t, i, j, k, p)
  ξ̂y(t, i, j, k, p) = z_η_x_ζ(t, i, j, k, p) − z_ζ_x_η(t, i, j, k, p)
  η̂y(t, i, j, k, p) = z_ζ_x_ξ(t, i, j, k, p) − z_ξ_x_ζ(t, i, j, k, p)
  ζ̂y(t, i, j, k, p) = z_ξ_x_η(t, i, j, k, p) − z_η_x_ξ(t, i, j, k, p)
  ξ̂z(t, i, j, k, p) = x_η_y_ζ(t, i, j, k, p) − x_ζ_y_η(t, i, j, k, p)
  η̂z(t, i, j, k, p) = x_ζ_y_ξ(t, i, j, k, p) − x_ξ_y_ζ(t, i, j, k, p)
  ζ̂z(t, i, j, k, p) = x_ξ_y_η(t, i, j, k, p) − x_η_y_ξ(t, i, j, k, p)

  ξ̂x_val_and_ξderivs = ξ_derivs(ξ̂x, backend, edge_scheme)
  η̂x_val_and_ξderivs = ξ_derivs(η̂x, backend, edge_scheme)
  ζ̂x_val_and_ξderivs = ξ_derivs(ζ̂x, backend, edge_scheme)
  ξ̂y_val_and_ξderivs = ξ_derivs(ξ̂y, backend, edge_scheme)
  η̂y_val_and_ξderivs = ξ_derivs(η̂y, backend, edge_scheme)
  ζ̂y_val_and_ξderivs = ξ_derivs(ζ̂y, backend, edge_scheme)
  ξ̂z_val_and_ξderivs = ξ_derivs(ξ̂z, backend, edge_scheme)
  η̂z_val_and_ξderivs = ξ_derivs(η̂z, backend, edge_scheme)
  ζ̂z_val_and_ξderivs = ξ_derivs(ζ̂z, backend, edge_scheme)

  ξ̂x_val_and_ηderivs = η_derivs(ξ̂x, backend, edge_scheme)
  η̂x_val_and_ηderivs = η_derivs(η̂x, backend, edge_scheme)
  ζ̂x_val_and_ηderivs = η_derivs(ζ̂x, backend, edge_scheme)
  ξ̂y_val_and_ηderivs = η_derivs(ξ̂y, backend, edge_scheme)
  η̂y_val_and_ηderivs = η_derivs(η̂y, backend, edge_scheme)
  ζ̂y_val_and_ηderivs = η_derivs(ζ̂y, backend, edge_scheme)
  ξ̂z_val_and_ηderivs = η_derivs(ξ̂z, backend, edge_scheme)
  η̂z_val_and_ηderivs = η_derivs(η̂z, backend, edge_scheme)
  ζ̂z_val_and_ηderivs = η_derivs(ζ̂z, backend, edge_scheme)

  ξ̂x_val_and_ζderivs = ζ_derivs(ξ̂x, backend, edge_scheme)
  η̂x_val_and_ζderivs = ζ_derivs(η̂x, backend, edge_scheme)
  ζ̂x_val_and_ζderivs = ζ_derivs(ζ̂x, backend, edge_scheme)
  ξ̂y_val_and_ζderivs = ζ_derivs(ξ̂y, backend, edge_scheme)
  η̂y_val_and_ζderivs = ζ_derivs(η̂y, backend, edge_scheme)
  ζ̂y_val_and_ζderivs = ζ_derivs(ζ̂y, backend, edge_scheme)
  ξ̂z_val_and_ζderivs = ζ_derivs(ξ̂z, backend, edge_scheme)
  η̂z_val_and_ζderivs = ζ_derivs(η̂z, backend, edge_scheme)
  ζ̂z_val_and_ζderivs = ζ_derivs(ζ̂z, backend, edge_scheme)

  Jinv_val_and_ξderivs = ξ_derivs(inverse_jacobian_matrix, backend, edge_scheme)
  Jinv_val_and_ηderivs = η_derivs(inverse_jacobian_matrix, backend, edge_scheme)
  Jinv_val_and_ζderivs = ζ_derivs(inverse_jacobian_matrix, backend, edge_scheme)

  ξ̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂x_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  η̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂x_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ζ̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂x_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ξ̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂y_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  η̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂y_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ζ̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂y_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ξ̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂z_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  η̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂z_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  ζ̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂z_val_and_ξderivs, t, i, j, k, p, edge_scheme)

  ξ̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂x_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  η̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂x_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ζ̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂x_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ξ̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂y_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  η̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂y_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ζ̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂y_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ξ̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂z_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  η̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂z_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  ζ̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂z_val_and_ηderivs, t, i, j, k, p, edge_scheme)

  ξ̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂x_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  η̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂x_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ζ̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂x_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ξ̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂y_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  η̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂y_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ζ̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂y_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ξ̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂z_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  η̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂z_val_and_ζderivs, t, i, j, k, p, edge_scheme)
  ζ̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂z_val_and_ζderivs, t, i, j, k, p, edge_scheme)

  Jinvᵢ₊½(t, i, j, k, p) = ϕ_iedge(Jinv_val_and_ξderivs, t, i, j, k, p, edge_scheme)
  Jinvⱼ₊½(t, i, j, k, p) = ϕ_jedge(Jinv_val_and_ηderivs, t, i, j, k, p, edge_scheme)
  Jinvₖ₊½(t, i, j, k, p) = ϕ_kedge(Jinv_val_and_ζderivs, t, i, j, k, p, edge_scheme)

  return (
    (;
      ξ̂x,
      η̂x,
      ζ̂x,
      ξ̂y,
      η̂y,
      ζ̂y,
      ξ̂z,
      η̂z,
      ζ̂z,
      # J=jacobian_matrix,
      Jinv=inverse_jacobian_matrix,
      # normJinv=normalized_inverse_jacobian_matrix,
    ),
    (;
      ξ̂xᵢ₊½,
      η̂xᵢ₊½,
      ζ̂xᵢ₊½,
      ξ̂yᵢ₊½,
      η̂yᵢ₊½,
      ζ̂yᵢ₊½,
      ξ̂zᵢ₊½,
      η̂zᵢ₊½,
      ζ̂zᵢ₊½,
      Jinvᵢ₊½,
      #
      ξ̂xⱼ₊½,
      η̂xⱼ₊½,
      ζ̂xⱼ₊½,
      ξ̂yⱼ₊½,
      η̂yⱼ₊½,
      ζ̂yⱼ₊½,
      ξ̂zⱼ₊½,
      η̂zⱼ₊½,
      ζ̂zⱼ₊½,
      Jinvⱼ₊½,
      #
      ξ̂xₖ₊½,
      η̂xₖ₊½,
      ζ̂xₖ₊½,
      ξ̂yₖ₊½,
      η̂yₖ₊½,
      ζ̂yₖ₊½,
      ξ̂zₖ₊½,
      η̂zₖ₊½,
      ζ̂zₖ₊½,
      Jinvₖ₊½,
    ),
  )
end
#-------------------------------------------------------------
#-------------------------------------------------------------

@inline _ξ_eval_3d_metriccache(ξ, ϕ, t, j, k, p) = ϕ(t, ξ, j, k, p)
@inline _η_eval_3d_metriccache(η, ϕ, t, i, k, p) = ϕ(t, i, η, k, p)
@inline _ζ_eval_3d_metriccache(ζ, ϕ, t, i, j, p) = ϕ(t, i, j, ζ, p)

ξ_derivs(ϕ, backend) = ξ_derivs(ϕ, backend, EdgeInterpolationOrder3())
η_derivs(ϕ, backend) = η_derivs(ϕ, backend, EdgeInterpolationOrder3())
ζ_derivs(ϕ, backend) = ζ_derivs(ϕ, backend, EdgeInterpolationOrder3())

function ξ_derivs(ϕ, backend, ::EdgeInterpolationOrder1)
  ϕval(t, i, j, k, p) = _ξ_eval_3d_metriccache(i, ϕ, t, j, k, p)
  return ϕval
end

function ξ_derivs(ϕ, backend, ::EdgeInterpolationOrder2)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_derivative(
    _ξ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_and_derivative(
      _ξ_eval_3d_metriccache,
      prep,
      backend,
      i,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function ξ_derivs(ϕ, backend, ::EdgeInterpolationOrder3)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_second_derivative(
    _ξ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(
      _ξ_eval_3d_metriccache,
      prep,
      backend,
      i,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function η_derivs(ϕ, backend, ::EdgeInterpolationOrder1)
  ϕval(t, i, j, k, p) = _η_eval_3d_metriccache(j, ϕ, t, i, k, p)
  return ϕval
end

function η_derivs(ϕ, backend, ::EdgeInterpolationOrder2)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_derivative(
    _η_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_and_derivative(
      _η_eval_3d_metriccache,
      prep,
      backend,
      j,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function η_derivs(ϕ, backend, ::EdgeInterpolationOrder3)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_second_derivative(
    _η_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(
      _η_eval_3d_metriccache,
      prep,
      backend,
      j,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(k),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function ζ_derivs(ϕ, backend, ::EdgeInterpolationOrder1)
  ϕval(t, i, j, k, p) = _ζ_eval_3d_metriccache(k, ϕ, t, i, j, p)
  return ϕval
end

function ζ_derivs(ϕ, backend, ::EdgeInterpolationOrder2)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_derivative(
    _ζ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_and_derivative(
      _ζ_eval_3d_metriccache,
      prep,
      backend,
      k,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end

function ζ_derivs(ϕ, backend, ::EdgeInterpolationOrder3)
  cϕ = DifferentiationInterface.Constant(ϕ)
  prep = prepare_second_derivative(
    _ζ_eval_3d_metriccache,
    backend,
    0.0,
    cϕ,
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(0.0),
    DifferentiationInterface.Constant(nothing);
    strict=Val(false),
  )

  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(
      _ζ_eval_3d_metriccache,
      prep,
      backend,
      k,
      cϕ,
      DifferentiationInterface.Constant(t),
      DifferentiationInterface.Constant(i),
      DifferentiationInterface.Constant(j),
      DifferentiationInterface.Constant(p),
    )
  end

  return ϕall
end
