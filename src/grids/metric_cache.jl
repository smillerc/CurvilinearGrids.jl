struct MetricCache{F1,F2}
  forward::F1
  inverse::F2
end

"""
3D metric cache
"""
function MetricCache(x::Function, y::Function, z::Function, backend)
  xξ(i, j, k) = derivative(ξ -> x(ξ, j, k), backend, i)
  xη(i, j, k) = derivative(η -> x(i, η, k), backend, j)
  xζ(i, j, k) = derivative(ζ -> x(i, j, ζ), backend, k)
  yξ(i, j, k) = derivative(ξ -> y(ξ, j, k), backend, i)
  yη(i, j, k) = derivative(η -> y(i, η, k), backend, j)
  yζ(i, j, k) = derivative(ζ -> y(i, j, ζ), backend, k)
  zξ(i, j, k) = derivative(ξ -> z(ξ, j, k), backend, i)
  zη(i, j, k) = derivative(η -> z(i, η, k), backend, j)
  zζ(i, j, k) = derivative(ζ -> z(i, j, ζ), backend, k)

  function jacobian_matrix(i, j, k)

    # compute the jacobian matrix w/o any extra logic
    function compute_jacobian_matrix(i, j, k)
      jac = DifferentiationInterface.jacobian(
        u -> SVector(x(u...), y(u...), z(u...)), backend, @SVector [i, j, k]
      )

      return jac
    end

    jac_matrix = compute_jacobian_matrix(i, j, k)
    J = det(jac_matrix)

    # sometimes the Jacobian can be zero for various reasons 
    # (derivs are zero)

    if !iszero(J) && isfinite(J)
      return jac_matrix
    else

      # perturb the coordinates if we're at a singularity

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
        ξηζ = (i, j, k) .+ pert

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

  jacobian(i, j, k) = det(jacobian_matrix(i, j, k))
  function jinv(i, j, k)
    J_inv = inv(jacobian_matrix(i, j, k))

    # J_inv = @. J_inv * (abs(J_inv) >= eps())
    return J_inv
  end

  function normalized_jinv(i, j, k)
    jac = jacobian_matrix(i, j, k)
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
    yξ=yξ,
    yη=yη,
    yζ=yζ,
    zξ=zξ,
    zη=zη,
    zζ=zζ,
  )

  x_ξ_y(i, j, k) = xξ(i, j, k) * y(i, j, k)
  x_η_y(i, j, k) = xη(i, j, k) * y(i, j, k)
  x_ζ_y(i, j, k) = xζ(i, j, k) * y(i, j, k)

  y_ξ_z(i, j, k) = yξ(i, j, k) * z(i, j, k)
  y_η_z(i, j, k) = yη(i, j, k) * z(i, j, k)
  y_ζ_z(i, j, k) = yζ(i, j, k) * z(i, j, k)

  z_ξ_x(i, j, k) = zξ(i, j, k) * x(i, j, k)
  z_η_x(i, j, k) = zη(i, j, k) * x(i, j, k)
  z_ζ_x(i, j, k) = zζ(i, j, k) * x(i, j, k)

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
  ξ̂x(i, j, k) = y_η_z_ζ(i, j, k) − y_ζ_z_η(i, j, k)
  η̂x(i, j, k) = y_ζ_z_ξ(i, j, k) − y_ξ_z_ζ(i, j, k)
  ζ̂x(i, j, k) = y_ξ_z_η(i, j, k) − y_η_z_ξ(i, j, k)
  ξ̂y(i, j, k) = z_η_x_ζ(i, j, k) − z_ζ_x_η(i, j, k)
  η̂y(i, j, k) = z_ζ_x_ξ(i, j, k) − z_ξ_x_ζ(i, j, k)
  ζ̂y(i, j, k) = z_ξ_x_η(i, j, k) − z_η_x_ξ(i, j, k)
  ξ̂z(i, j, k) = x_η_y_ζ(i, j, k) − x_ζ_y_η(i, j, k)
  η̂z(i, j, k) = x_ζ_y_ξ(i, j, k) − x_ξ_y_ζ(i, j, k)
  ζ̂z(i, j, k) = x_ξ_y_η(i, j, k) − x_η_y_ξ(i, j, k)

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
  xξ(i, j) = derivative(ξ -> x(ξ, j), backend, i)
  xη(i, j) = derivative(η -> x(i, η), backend, j)
  yξ(i, j) = derivative(ξ -> y(ξ, j), backend, i)
  yη(i, j) = derivative(η -> y(i, η), backend, j)

  function jacobian_matrix(i, j)

    # compute the jacobian matrix w/o any extra logic
    function compute_jacobian_matrix(i, j)
      jac = DifferentiationInterface.jacobian(
        u -> SVector(x(u...), y(u...)), backend, @SVector [i, j]
      )

      return jac
    end

    jac_matrix = compute_jacobian_matrix(i, j)
    J = det(jac_matrix)

    # sometimes the Jacobian can be zero for various reasons 
    # (derivs are zero)

    if !iszero(J) && isfinite(J)
      return jac_matrix
    else
      # perturb the coordinates if we're at a singularity

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
        ξηζ = (i, j) .+ pert

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

  jacobian(i, j) = det(jacobian_matrix(i, j))
  function jinv(i, j)
    J_inv = inv(jacobian_matrix(i, j))
    return J_inv
  end

  function normalized_jinv(i, j)
    jac = jacobian_matrix(i, j)
    J = det(jac)
    J_inv = inv(jac)

    return J_inv .* J
  end

  forward_metrics = (; jacobian=jacobian_matrix, J=jacobian, xξ=xξ, xη=xη, yξ=yξ, yη=yη)

  #   x_ξ_y(i, j) = xξ(i, j) * y(i, j)
  #   x_η_y(i, j) = xη(i, j) * y(i, j)

  #   y_ξ_x(i, j) = yξ(i, j) * x(i, j)
  #   y_η_x(i, j) = yη(i, j) * x(i, j)

  #   _, x_ξ_y_η = cell_center_derivative_2d(x_ξ_y, backend)
  #   #   x_η_y_ξ, _, = cell_center_derivative_2d(x_η_y, backend)

  #   _, y_ξ_x_η = cell_center_derivative_2d(y_ξ_x, backend)
  #   #   y_η_x_ξ, _, = cell_center_derivative_2d(y_η_x, backend)

  #   ξ̂x(i, j) = yη(i, j) / (x_ξ_y_η(i, j) − y_ξ_x_η(i, j))
  #   ξ̂y(i, j) = -xη(i, j) / (x_ξ_y_η(i, j) − y_ξ_x_η(i, j))
  #   η̂x(i, j) = -yξ(i, j) / (x_ξ_y_η(i, j) − y_ξ_x_η(i, j))
  #   η̂y(i, j) = xξ(i, j) / (x_ξ_y_η(i, j) − y_ξ_x_η(i, j))

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
  xξ(i) = derivative(ξ -> x(ξ), backend, i)

  jacobian_matrix(i) = @SMatrix [xξ(i)]

  jacobian(i) = det(jacobian_matrix(i))
  jinv(i) = inv(jacobian_matrix(i))

  function normalized_jinv(i)
    jac = jacobian_matrix(i)
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

  ∂ϕ_∂ξ(i, j, k) = ϕᵢ₊½(i, j, k) - ϕᵢ₊½(i - 1, j, k)
  ∂ϕ_∂η(i, j, k) = ϕⱼ₊½(i, j, k) - ϕⱼ₊½(i, j - 1, k)
  ∂ϕ_∂ζ(i, j, k) = ϕₖ₊½(i, j, k) - ϕₖ₊½(i, j, k - 1)

  return (; ∂ϕ_∂ξ, ∂ϕ_∂η, ∂ϕ_∂ζ)
end

function cell_center_derivative_2d(ϕ, backend)
  ϕᵢ₊½, ϕⱼ₊½ = edge_functions_2d(ϕ, backend)

  ∂ϕ_∂ξ(i, j) = ϕᵢ₊½(i, j) - ϕᵢ₊½(i - 1, j)
  ∂ϕ_∂η(i, j) = ϕⱼ₊½(i, j) - ϕⱼ₊½(i, j - 1)

  return (; ∂ϕ_∂ξ, ∂ϕ_∂η)
end

function cell_center_derivative_1d(ϕ, backend)
  ϕᵢ₊½ = edge_functions_1d(ϕ, backend)

  ∂ϕ_∂ξ(i) = ϕᵢ₊½(i) - ϕᵢ₊½(i - 1)

  return (; ∂ϕ_∂ξ,)
end

function edge_functions_3d(ϕ, backend)

  # returns val, ∂, ∂²
  ξ_derivs(i, j, k) = value_derivative_and_second_derivative(ξ -> ϕ(ξ, j, k), backend, i)
  η_derivs(i, j, k) = value_derivative_and_second_derivative(η -> ϕ(i, η, k), backend, j)
  ζ_derivs(i, j, k) = value_derivative_and_second_derivative(ζ -> ϕ(i, j, ζ), backend, k)

  function ϕᵢ₊½(i, j, k)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(i, j, k)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(i + 1, j, k)

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ϕⱼ₊½(i, j, k)
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = η_derivs(i, j, k)
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = η_derivs(i, j + 1, k)

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  function ϕₖ₊½(i, j, k)
    ϕₖ, ∂ϕ_∂ξₖ, ∂²ϕ_∂ξ²ₖ = ζ_derivs(i, j, k)
    ϕₖ₊₁, ∂ϕ_∂ξₖ₊₁, ∂²ϕ_∂ξ²ₖ₊₁ = ζ_derivs(i, j, k + 1)

    ϕᴸₖ₊½ = ϕₖ + (1 / 2) * ∂ϕ_∂ξₖ + (1 / 12) * ∂²ϕ_∂ξ²ₖ
    ϕᴿₖ₊½ = ϕₖ₊₁ - (1 / 2) * ∂ϕ_∂ξₖ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ₖ₊₁

    return (ϕᴸₖ₊½ + ϕᴿₖ₊½) / 2
  end

  return (; ϕᵢ₊½, ϕⱼ₊½, ϕₖ₊½)
end

function edge_functions_2d(ϕ, backend)

  # returns val, ∂, ∂²
  ξ_derivs(i, j) = value_derivative_and_second_derivative(ξ -> ϕ(ξ, j), backend, i)
  η_derivs(i, j) = value_derivative_and_second_derivative(η -> ϕ(i, η), backend, j)

  function ϕᵢ₊½(i, j)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(i, j)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(i + 1, j)

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  function ϕⱼ₊½(i, j)
    ϕⱼ, ∂ϕ_∂ξⱼ, ∂²ϕ_∂ξ²ⱼ = η_derivs(i, j)
    ϕⱼ₊₁, ∂ϕ_∂ξⱼ₊₁, ∂²ϕ_∂ξ²ⱼ₊₁ = η_derivs(i, j + 1)

    ϕᴸⱼ₊½ = ϕⱼ + (1 / 2) * ∂ϕ_∂ξⱼ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ
    ϕᴿⱼ₊½ = ϕⱼ₊₁ - (1 / 2) * ∂ϕ_∂ξⱼ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ⱼ₊₁

    return (ϕᴸⱼ₊½ + ϕᴿⱼ₊½) / 2
  end

  return (; ϕᵢ₊½, ϕⱼ₊½)
end

function edge_functions_1d(ϕ, backend)

  # returns val, ∂, ∂²
  ξ_derivs(i) = value_derivative_and_second_derivative(ξ -> ϕ(ξ), backend, i)

  function ϕᵢ₊½(i)
    ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ξ_derivs(i)
    ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ξ_derivs(i + 1)

    ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
    ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

    return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
  end

  return (; ϕᵢ₊½)
end
