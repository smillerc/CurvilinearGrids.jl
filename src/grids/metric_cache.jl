struct MetricCache{F1,F2}
  forward::F1
  inverse::F2
end

"""
3D metric cache
"""
function MetricCache(x::Function, y::Function, z::Function, backend)
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
    J = det(jac_matrix)

    # sometimes the Jacobian can be zero for various reasons 
    # (derivs are zero)

    if !iszero(J) && isfinite(J)
      return jac_matrix
    else

      # perturb the coordinates if we're at a singularity
      @warn "perturbing the coorindates to avoid a Jacobian singularity"
      pert = sqrt(eps())
      iter = 1
      itermax = 25

      last_valid_jacobian_matrix = jac_matrix

      # perturb until we get zero again, and then use the last nonzero value
      while true
        if iter > itermax
          error(
            "Maximum iteration count reached ($itermax) when trying to iteratively calculate the Jacobian near a singularity",
          )
        end

        # perturb the coordinates
        ξηζ = (t, i .+ pert, j .+ pert, k .+ pert, p)

        # compute the jacobian
        _jac = compute_jacobian_matrix(ξηζ...)
        _J = det(_jac)

        if !isfinite(_J) || iszero(_J)
          # we've made the perturbation too small
          break
        else
          last_valid_jacobian_matrix = _jac
        end

        # make the perturbation smaller
        pert = pert / 10
        iter += 1
      end

      last_valid_jacobian_matrix = @. last_valid_jacobian_matrix *
        (abs(last_valid_jacobian_matrix) >= eps())
      return last_valid_jacobian_matrix
    end
  end

  jacobian(t, i, j, k, p) = det(jacobian_matrix(t, i, j, k, p))
  function jinv(t, i, j, k, p)
    J_inv = inv(jacobian_matrix(t, i, j, k, p))

    # J_inv = @. J_inv * (abs(J_inv) >= eps())
    return J_inv
  end

  function normalized_jinv(t, i, j, k, p)
    jac = jacobian_matrix(t, i, j, k, p)
    J = det(jac)
    J_inv = inv(jac)

    return J_inv .* J
  end

  forward_metrics = (;
    jacobian=jacobian_matrix,
    J=jacobian,
    xξ=xξ,
    xη=xη,
    xζ=xζ,
    xτ=xτ,
    yξ=yξ,
    yη=yη,
    yζ=yζ,
    yτ=yτ,
    zξ=zξ,
    zη=zη,
    zζ=zζ,
    zτ=zτ,
  )

  x_ξ_y(t, i, j, k, p) = xξ(t, i, j, k, p) * y(t, i, j, k, p)
  x_η_y(t, i, j, k, p) = xη(t, i, j, k, p) * y(t, i, j, k, p)
  x_ζ_y(t, i, j, k, p) = xζ(t, i, j, k, p) * y(t, i, j, k, p)

  y_ξ_z(t, i, j, k, p) = yξ(t, i, j, k, p) * z(t, i, j, k, p)
  y_η_z(t, i, j, k, p) = yη(t, i, j, k, p) * z(t, i, j, k, p)
  y_ζ_z(t, i, j, k, p) = yζ(t, i, j, k, p) * z(t, i, j, k, p)

  z_ξ_x(t, i, j, k, p) = zξ(t, i, j, k, p) * x(t, i, j, k, p)
  z_η_x(t, i, j, k, p) = zη(t, i, j, k, p) * x(t, i, j, k, p)
  z_ζ_x(t, i, j, k, p) = zζ(t, i, j, k, p) * x(t, i, j, k, p)

  _, x_ξ_y_η, x_ξ_y_ζ = cell_center_derivative_3d(x_ξ_y, backend)
  x_η_y_ξ, _, x_η_y_ζ = cell_center_derivative_3d(x_η_y, backend)
  x_ζ_y_ξ, x_ζ_y_η, _ = cell_center_derivative_3d(x_ζ_y, backend)

  _, y_ξ_z_η, y_ξ_z_ζ = cell_center_derivative_3d(y_ξ_z, backend)
  y_η_z_ξ, _, y_η_z_ζ = cell_center_derivative_3d(y_η_z, backend)
  y_ζ_z_ξ, y_ζ_z_η, _ = cell_center_derivative_3d(y_ζ_z, backend)

  _, z_ξ_x_η, z_ξ_x_ζ = cell_center_derivative_3d(z_ξ_x, backend)
  z_η_x_ξ, _, z_η_x_ζ = cell_center_derivative_3d(z_η_x, backend)
  z_ζ_x_ξ, z_ζ_x_η, _ = cell_center_derivative_3d(z_ζ_x, backend)

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

  inverse_metrics = (;
    ξ̂x=ξ̂x,
    ξ̂y=ξ̂y,
    ξ̂z=ξ̂z,
    η̂x=η̂x,
    η̂y=η̂y,
    η̂z=η̂z,
    ζ̂x=ζ̂x,
    ζ̂y=ζ̂y,
    ζ̂z=ζ̂z,
    Jinv=jinv,
    Jinv_norm=normalized_jinv,
  )

  return MetricCache(forward_metrics, inverse_metrics)
end

"""
2D metric cache
"""
function MetricCache(x::Function, y::Function, backend)
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
    # end

    # jac_matrix = compute_jacobian_matrix(t, i, j, p)
    # J = det(jac_matrix)

    # # sometimes the Jacobian can be zero for various reasons 
    # # (derivs are zero)

    # if !iszero(J) && isfinite(J)
    #   return jac_matrix
    # else
    #   # perturb the coordinates if we're at a singularity

    #   pert = sqrt(eps())
    #   iter = 1
    #   itermax = 25

    #   last_valid_jacobian_matrix = jac_matrix

    #   # perturb until we get zero again, and then use the last nonzero value
    #   while true
    #     if iter > itermax
    #       error(
    #         "Maximum iteration count reached ($itermax) when trying to iteratively calculate the Jacobian near a singularity",
    #       )
    #     end

    #     # perturb the coordinates
    #     ξη = (i, j) .+ pert

    #     # compute the jacobian
    #     _jac = compute_jacobian_matrix(t, ξη..., p)
    #     _J = det(_jac)

    #     if !isfinite(_J) || iszero(_J)
    #       # we've made the perturbation too small
    #       break
    #     else
    #       last_valid_jacobian_matrix = _jac
    #     end

    #     # make the perturbation smaller
    #     pert = pert / 10
    #     iter += 1
    #   end

    #   last_valid_jacobian_matrix = @. last_valid_jacobian_matrix *
    #     (abs(last_valid_jacobian_matrix) >= eps())
    #   return last_valid_jacobian_matrix
    # end
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

  #   x_ξ_y(t, i, j, p) = xξ(t, i, j, p) * y(t, i, j, p)
  #   x_η_y(t, i, j, p) = xη(t, i, j, p) * y(t, i, j, p)

  #   y_ξ_x(t, i, j, p) = yξ(t, i, j, p) * x(t, i, j, p)
  #   y_η_x(t, i, j, p) = yη(t, i, j, p) * x(t, i, j, p)

  #   _, x_ξ_y_η = cell_center_derivative_2d(x_ξ_y, backend)
  #   #   x_η_y_ξ, _, = cell_center_derivative_2d(x_η_y, backend)

  #   _, y_ξ_x_η = cell_center_derivative_2d(y_ξ_x, backend)
  #   #   y_η_x_ξ, _, = cell_center_derivative_2d(y_η_x, backend)

  #   ξ̂x(t, i, j, p) = yη(t, i, j, p) / (x_ξ_y_η(t, i, j, p) − y_ξ_x_η(t, i, j, p))
  #   ξ̂y(t, i, j, p) = -xη(t, i, j, p) / (x_ξ_y_η(t, i, j, p) − y_ξ_x_η(t, i, j, p))
  #   η̂x(t, i, j, p) = -yξ(t, i, j, p) / (x_ξ_y_η(t, i, j, p) − y_ξ_x_η(t, i, j, p))
  #   η̂y(t, i, j, p) = xξ(t, i, j, p) / (x_ξ_y_η(t, i, j, p) − y_ξ_x_η(t, i, j, p))

  inverse_metrics = (;
    # ξ̂x=ξ̂x, ξ̂y=ξ̂y, η̂x=η̂x, η̂y=η̂y,
    Jinv=jinv,
    Jinv_norm=normalized_jinv,
  )

  return MetricCache(forward_metrics, inverse_metrics)
end

"""
1D metric cache
"""
function MetricCache(x::Function, backend)
  xξ(t, i, p) = derivative(ξ -> x(t, ξ, p), backend, i)
  xτ(t, i, p) = derivative(τ -> x(τ, i, p), backend, t)

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
  return MetricCache(forward_metrics, inverse_metrics)
end

function cell_center_derivative_3d(ϕ, backend)
  ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½ = edge_functions_3d(ϕ, backend)

  ∂ϕ_∂ξ(t, i, j, k, p) = ϕᵢ₊½(t, i, j, k, p) - ϕᵢ₊½(t, i - 1, j, k, p)
  ∂ϕ_∂η(t, i, j, k, p) = ϕⱼ₊½(t, i, j, k, p) - ϕⱼ₊½(t, i, j - 1, k, p)
  ∂ϕ_∂ζ(t, i, j, k, p) = ϕₖ₊½(t, i, j, k, p) - ϕₖ₊½(t, i, j, k - 1, p)

  return (; ∂ϕ_∂ξ, ∂ϕ_∂η, ∂ϕ_∂ζ)
end

function cell_center_derivative_2d(ϕ, backend)
  ϕᵢ₊½, ϕⱼ₊½ = edge_functions_2d(ϕ, backend)

  ∂ϕ_∂ξ(t, i, j, p) = ϕᵢ₊½(t, i, j, p) - ϕᵢ₊½(t, i - 1, j, p)
  ∂ϕ_∂η(t, i, j, p) = ϕⱼ₊½(t, i, j, p) - ϕⱼ₊½(t, i, j - 1, p)

  return (; ∂ϕ_∂ξ, ∂ϕ_∂η)
end

function cell_center_derivative_1d(ϕ, backend)
  ϕᵢ₊½ = edge_functions_1d(ϕ, backend)

  ∂ϕ_∂ξ(t, i, p) = ϕᵢ₊½(t, i, p) - ϕᵢ₊½(t, i - 1, p)

  return (; ∂ϕ_∂ξ,)
end

function edge_functions_3d(ϕ, backend)

  # returns val, ∂, ∂²
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

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ϕⱼ₊½(t, i, j, k, p)
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = η_derivs(t, i, j, k, p)
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = η_derivs(t, i, j + 1, k, p)

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  function ϕₖ₊½(t, i, j, k, p)
    ϕₖ, ∂ϕ_∂ξₖ, ∂²ϕ_∂ξ²ₖ = ζ_derivs(t, i, j, k, p)
    ϕₖ₊₁, ∂ϕ_∂ξₖ₊₁, ∂²ϕ_∂ξ²ₖ₊₁ = ζ_derivs(t, i, j, k + 1, p)

    ϕᴸₖ₊½ = ϕₖ + (1 / 2) * ∂ϕ_∂ξₖ + (1 / 12) * ∂²ϕ_∂ξ²ₖ
    ϕᴿₖ₊½ = ϕₖ₊₁ - (1 / 2) * ∂ϕ_∂ξₖ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ₖ₊₁

    return (ϕᴸₖ₊½ + ϕᴿₖ₊½) / 2
  end

  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

function edge_functions_2d(ϕ, backend)

  # returns val, ∂, ∂²
  function ξ_derivs(t, i, j, p)
    value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, j, p), backend, i)
  end
  function η_derivs(t, i, j, p)
    value_derivative_and_second_derivative(η -> ϕ(t, i, η, p), backend, j)
  end

  function ϕᵢ₊½(t, i, j, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, j, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, j, p)

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ϕⱼ₊½(t, i, j, p)
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = η_derivs(t, i, j, p)
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = η_derivs(t, i, j + 1, p)

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  return (; ϕᵢ₊½, ϕⱼ₊½)
end

function edge_functions_1d(ϕ, backend)

  # returns val, ∂, ∂²
  ξ_derivs(t, i, p) = value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, p), backend, i)

  function ϕᵢ₊½(t, i, p)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(t, i, p)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(t, i + 1, p)

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  return (; ϕᵢ₊½)
end
