using DifferentiationInterface
using ForwardDiff: ForwardDiff
using UnPack
using CartesianDomains
using BenchmarkTools
using Polyester

function wavy_mapping()
  function x(t, i, j, k, p)
    @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin = p
    xmin + Δx0 * ((i - 1) + Ax * sin(pi * n * (j - 1) * Δy0) * sin(pi * n * (k - 1) * Δz0))
  end

  function y(t, i, j, k, p)
    @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, ymin = p
    ymin + Δy0 * ((j - 1) + Ay * sin(pi * n * (k - 1) * Δz0) * sin(pi * n * (i - 1) * Δx0))
  end

  function z(t, i, j, k, p)
    @unpack Ax, Ay, Az, n, Δx0, Δy0, Δz0, zmin = p
    zmin + Δz0 * ((k - 1) + Az * sin(pi * n * (i - 1) * Δx0) * sin(pi * n * (j - 1) * Δy0))
  end

  return (x, y, z)
end

function wavy_params(ncells)
  ni, nj, nk = ncells
  Lx = Ly = Lz = 12

  xmin = -Lx / 2
  ymin = -Ly / 2
  zmin = -Lz / 2

  Δx0 = Lx / ni
  Δy0 = Ly / nj
  Δz0 = Lz / nk

  Ax = 0.2 / Δx0
  Ay = 0.4 / Δy0
  Az = 0.6 / Δz0

  n = 0.5

  params = (; Ax, Ay, Az, n, Δx0, Δy0, Δz0, xmin, ymin, zmin)
end

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

function get_ξ̂x(x, y, z, backend)
  yη(t, i, j, k, p) = derivative(η -> y(t, i, η, k, p), backend, j)
  yζ(t, i, j, k, p) = derivative(ζ -> y(t, i, j, ζ, p), backend, k)

  y_η_z(t, i, j, k, p) = yη(t, i, j, k, p) * z(t, i, j, k, p)
  y_ζ_z(t, i, j, k, p) = yζ(t, i, j, k, p) * z(t, i, j, k, p)

  y_η_z_ζ = ∂ϕ_∂ζ_3d(y_η_z, backend)
  y_ζ_z_η = ∂ϕ_∂η_3d(y_ζ_z, backend)

  ξ̂x(t, i, j, k, p) = y_η_z_ζ(t, i, j, k, p) − y_ζ_z_η(t, i, j, k, p)

  return ξ̂x
end

function get_η̂x(x, y, z, backend)
  yζ(t, i, j, k, p) = derivative(ζ -> y(t, i, j, ζ, p), backend, k)
  yξ(t, i, j, k, p) = derivative(ξ -> y(t, ξ, j, k, p), backend, i)

  y_ζ_z(t, i, j, k, p) = yζ(t, i, j, k, p) * z(t, i, j, k, p)
  y_ξ_z(t, i, j, k, p) = yξ(t, i, j, k, p) * z(t, i, j, k, p)

  y_ζ_z_ξ = ∂ϕ_∂ξ_3d(y_ζ_z, backend)
  y_ξ_z_ζ = ∂ϕ_∂ζ_3d(y_ξ_z, backend)

  η̂x(t, i, j, k, p) = y_ζ_z_ξ(t, i, j, k, p) − y_ξ_z_ζ(t, i, j, k, p)
  return η̂x
end

function get_ζ̂x(x, y, z, backend)
  yξ(t, i, j, k, p) = derivative(ξ -> y(t, ξ, j, k, p), backend, i)
  yη(t, i, j, k, p) = derivative(η -> y(t, i, η, k, p), backend, j)

  y_ξ_z(t, i, j, k, p) = yξ(t, i, j, k, p) * z(t, i, j, k, p)
  y_η_z(t, i, j, k, p) = yη(t, i, j, k, p) * z(t, i, j, k, p)

  y_ξ_z_η = ∂ϕ_∂η_3d(y_ξ_z, backend)
  y_η_z_ξ = ∂ϕ_∂ξ_3d(y_η_z, backend)

  ζ̂x(t, i, j, k, p) = y_ξ_z_η(t, i, j, k, p) − y_η_z_ξ(t, i, j, k, p)
  return ζ̂x
end

function get_ξ̂y(x, y, z, backend)
  zη(t, i, j, k, p) = derivative(η -> z(t, i, η, k, p), backend, j)
  zζ(t, i, j, k, p) = derivative(ζ -> z(t, i, j, ζ, p), backend, k)

  z_η_x(t, i, j, k, p) = zη(t, i, j, k, p) * x(t, i, j, k, p)
  z_ζ_x(t, i, j, k, p) = zζ(t, i, j, k, p) * x(t, i, j, k, p)

  z_η_x_ζ = ∂ϕ_∂ζ_3d(z_η_x, backend)
  z_ζ_x_η = ∂ϕ_∂η_3d(z_ζ_x, backend)

  ξ̂y(t, i, j, k, p) = z_η_x_ζ(t, i, j, k, p) − z_ζ_x_η(t, i, j, k, p)
  return ξ̂y
end

function get_η̂y(x, y, z, backend)
  zζ(t, i, j, k, p) = derivative(ζ -> z(t, i, j, ζ, p), backend, k)
  zξ(t, i, j, k, p) = derivative(ξ -> z(t, ξ, j, k, p), backend, i)

  z_ζ_x(t, i, j, k, p) = zζ(t, i, j, k, p) * x(t, i, j, k, p)
  z_ξ_x(t, i, j, k, p) = zξ(t, i, j, k, p) * x(t, i, j, k, p)

  z_ζ_x_ξ = ∂ϕ_∂ξ_3d(z_ζ_x, backend)
  z_ξ_x_ζ = ∂ϕ_∂ζ_3d(z_ξ_x, backend)

  η̂y(t, i, j, k, p) = z_ζ_x_ξ(t, i, j, k, p) − z_ξ_x_ζ(t, i, j, k, p)
  return η̂y
end

function get_ζ̂y(x, y, z, backend)
  zξ(t, i, j, k, p) = derivative(ξ -> z(t, ξ, j, k, p), backend, i)
  zη(t, i, j, k, p) = derivative(η -> z(t, i, η, k, p), backend, j)

  z_ξ_x(t, i, j, k, p) = zξ(t, i, j, k, p) * x(t, i, j, k, p)
  z_η_x(t, i, j, k, p) = zη(t, i, j, k, p) * x(t, i, j, k, p)

  z_ξ_x_η = ∂ϕ_∂η_3d(z_ξ_x, backend)
  z_η_x_ξ = ∂ϕ_∂ξ_3d(z_η_x, backend)

  ζ̂y(t, i, j, k, p) = z_ξ_x_η(t, i, j, k, p) − z_η_x_ξ(t, i, j, k, p)
  return ζ̂y
end

function get_ξ̂z(x, y, z, backend)
  xη(t, i, j, k, p) = derivative(η -> x(t, i, η, k, p), backend, j)
  xζ(t, i, j, k, p) = derivative(ζ -> x(t, i, j, ζ, p), backend, k)

  x_η_y(t, i, j, k, p) = xη(t, i, j, k, p) * y(t, i, j, k, p)
  x_ζ_y(t, i, j, k, p) = xζ(t, i, j, k, p) * y(t, i, j, k, p)

  x_η_y_ζ = ∂ϕ_∂ζ_3d(x_η_y, backend)
  x_ζ_y_η = ∂ϕ_∂η_3d(x_ζ_y, backend)

  ξ̂z(t, i, j, k, p) = x_η_y_ζ(t, i, j, k, p) − x_ζ_y_η(t, i, j, k, p)
  return ξ̂z
end

function get_η̂z(x, y, z, backend)
  xζ(t, i, j, k, p) = derivative(ζ -> x(t, i, j, ζ, p), backend, k)
  xξ(t, i, j, k, p) = derivative(ξ -> x(t, i, j, k, p), backend, i)

  x_ζ_y(t, i, j, k, p) = xζ(t, i, j, k, p) * y(t, i, j, k, p)
  x_ξ_y(t, i, j, k, p) = xξ(t, i, j, k, p) * y(t, i, j, k, p)

  x_ζ_y_ξ = ∂ϕ_∂ξ_3d(x_ζ_y, backend)
  x_ξ_y_ζ = ∂ϕ_∂ζ_3d(x_ξ_y, backend)

  η̂z(t, i, j, k, p) = x_ζ_y_ξ(t, i, j, k, p) − x_ξ_y_ζ(t, i, j, k, p)
  return η̂z
end

function get_ζ̂z(x, y, z, backend)
  xξ(t, i, j, k, p) = derivative(ξ -> x(t, ξ, j, k, p), backend, i)
  xη(t, i, j, k, p) = derivative(η -> x(t, i, η, k, p), backend, j)

  x_ξ_y(t, i, j, k, p) = xξ(t, i, j, k, p) * y(t, i, j, k, p)
  x_η_y(t, i, j, k, p) = xη(t, i, j, k, p) * y(t, i, j, k, p)

  x_ξ_y_η = ∂ϕ_∂η_3d(x_ξ_y, backend)
  x_η_y_ξ = ∂ϕ_∂ξ_3d(x_η_y, backend)

  ζ̂z(t, i, j, k, p) = x_ξ_y_η(t, i, j, k, p) − x_η_y_ξ(t, i, j, k, p)
  return ζ̂z
end

function ξ_derivs(ϕ, backend)
  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(ξ -> ϕ(t, ξ, j, k, p), backend, i)
  end
  return ϕall
end

function η_derivs(ϕ, backend)
  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(η -> ϕ(t, i, η, k, p), backend, j)
  end
  return ϕall
end

function ζ_derivs(ϕ, backend)
  function ϕall(t, i, j, k, p)
    value_derivative_and_second_derivative(ζ -> ϕ(t, i, j, ζ, p), backend, k)
  end
  return ϕall
end

function ϕ_ᵢ₊½_3d(ϕ_set, t, i, j, k, p)
  ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ϕ_set(t, i, j, k, p)
  ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ϕ_set(t, i + 1, j, k, p)

  ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

  return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
end

function ϕ_iedge(ϕ_val_and_derivs, t, i, j, k, p)
  ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ϕ_val_and_derivs(t, i + 1, j, k, p)

  ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

  return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
end

function ϕ_jedge(ϕ_val_and_derivs, t, i, j, k, p)
  ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ϕ_val_and_derivs(t, i, j + 1, k, p)

  ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

  return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
end

function ϕ_kedge(ϕ_val_and_derivs, t, i, j, k, p)
  ϕᵢ, ∂ϕ_∂ξᵢ, ∂²ϕ_∂ξ²ᵢ = ϕ_val_and_derivs(t, i, j, k, p)
  ϕᵢ₊₁, ∂ϕ_∂ξᵢ₊₁, ∂²ϕ_∂ξ²ᵢ₊₁ = ϕ_val_and_derivs(t, i, j, k + 1, p)

  ϕᴸᵢ₊½ = ϕᵢ + (1 / 2) * ∂ϕ_∂ξᵢ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ
  ϕᴿᵢ₊½ = ϕᵢ₊₁ - (1 / 2) * ∂ϕ_∂ξᵢ₊₁ + (1 / 12) * ∂²ϕ_∂ξ²ᵢ₊₁

  return (ϕᴸᵢ₊½ + ϕᴿᵢ₊½) / 2
end

# @benchmark runme($xmap, $ymap, $zmap, $params, $backend)

function get_inverse_metric_terms(xmap, ymap, zmap, backend)
  ξ̂x = get_ξ̂x(xmap, ymap, zmap, backend)
  η̂x = get_η̂x(xmap, ymap, zmap, backend)
  ζ̂x = get_ζ̂x(xmap, ymap, zmap, backend)
  ξ̂y = get_ξ̂y(xmap, ymap, zmap, backend)
  η̂y = get_η̂y(xmap, ymap, zmap, backend)
  ζ̂y = get_ζ̂y(xmap, ymap, zmap, backend)
  ξ̂z = get_ξ̂z(xmap, ymap, zmap, backend)
  η̂z = get_η̂z(xmap, ymap, zmap, backend)
  ζ̂z = get_ζ̂z(xmap, ymap, zmap, backend)

  ξ̂x_val_and_ξderivs = ξ_derivs(ξ̂x, backend)
  η̂x_val_and_ξderivs = ξ_derivs(η̂x, backend)
  ζ̂x_val_and_ξderivs = ξ_derivs(ζ̂x, backend)
  ξ̂y_val_and_ξderivs = ξ_derivs(ξ̂y, backend)
  η̂y_val_and_ξderivs = ξ_derivs(η̂y, backend)
  ζ̂y_val_and_ξderivs = ξ_derivs(ζ̂y, backend)
  ξ̂z_val_and_ξderivs = ξ_derivs(ξ̂z, backend)
  η̂z_val_and_ξderivs = ξ_derivs(η̂z, backend)
  ζ̂z_val_and_ξderivs = ξ_derivs(ζ̂z, backend)

  ξ̂x_val_and_ηderivs = η_derivs(ξ̂x, backend)
  η̂x_val_and_ηderivs = η_derivs(η̂x, backend)
  ζ̂x_val_and_ηderivs = η_derivs(ζ̂x, backend)
  ξ̂y_val_and_ηderivs = η_derivs(ξ̂y, backend)
  η̂y_val_and_ηderivs = η_derivs(η̂y, backend)
  ζ̂y_val_and_ηderivs = η_derivs(ζ̂y, backend)
  ξ̂z_val_and_ηderivs = η_derivs(ξ̂z, backend)
  η̂z_val_and_ηderivs = η_derivs(η̂z, backend)
  ζ̂z_val_and_ηderivs = η_derivs(ζ̂z, backend)

  ξ̂x_val_and_ζderivs = ζ_derivs(ξ̂x, backend)
  η̂x_val_and_ζderivs = ζ_derivs(η̂x, backend)
  ζ̂x_val_and_ζderivs = ζ_derivs(ζ̂x, backend)
  ξ̂y_val_and_ζderivs = ζ_derivs(ξ̂y, backend)
  η̂y_val_and_ζderivs = ζ_derivs(η̂y, backend)
  ζ̂y_val_and_ζderivs = ζ_derivs(ζ̂y, backend)
  ξ̂z_val_and_ζderivs = ζ_derivs(ξ̂z, backend)
  η̂z_val_and_ζderivs = ζ_derivs(η̂z, backend)
  ζ̂z_val_and_ζderivs = ζ_derivs(ζ̂z, backend)

  ξ̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂x_val_and_ξderivs, t, i, j, k, p)
  η̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂x_val_and_ξderivs, t, i, j, k, p)
  ζ̂xᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂x_val_and_ξderivs, t, i, j, k, p)
  ξ̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂y_val_and_ξderivs, t, i, j, k, p)
  η̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂y_val_and_ξderivs, t, i, j, k, p)
  ζ̂yᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂y_val_and_ξderivs, t, i, j, k, p)
  ξ̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(ξ̂z_val_and_ξderivs, t, i, j, k, p)
  η̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(η̂z_val_and_ξderivs, t, i, j, k, p)
  ζ̂zᵢ₊½(t, i, j, k, p) = ϕ_iedge(ζ̂z_val_and_ξderivs, t, i, j, k, p)
  ξ̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂x_val_and_ηderivs, t, i, j, k, p)
  η̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂x_val_and_ηderivs, t, i, j, k, p)
  ζ̂xⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂x_val_and_ηderivs, t, i, j, k, p)
  ξ̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂y_val_and_ηderivs, t, i, j, k, p)
  η̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂y_val_and_ηderivs, t, i, j, k, p)
  ζ̂yⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂y_val_and_ηderivs, t, i, j, k, p)
  ξ̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(ξ̂z_val_and_ηderivs, t, i, j, k, p)
  η̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(η̂z_val_and_ηderivs, t, i, j, k, p)
  ζ̂zⱼ₊½(t, i, j, k, p) = ϕ_jedge(ζ̂z_val_and_ηderivs, t, i, j, k, p)
  ξ̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂x_val_and_ζderivs, t, i, j, k, p)
  η̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂x_val_and_ζderivs, t, i, j, k, p)
  ζ̂xₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂x_val_and_ζderivs, t, i, j, k, p)
  ξ̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂y_val_and_ζderivs, t, i, j, k, p)
  η̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂y_val_and_ζderivs, t, i, j, k, p)
  ζ̂yₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂y_val_and_ζderivs, t, i, j, k, p)
  ξ̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(ξ̂z_val_and_ζderivs, t, i, j, k, p)
  η̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(η̂z_val_and_ζderivs, t, i, j, k, p)
  ζ̂zₖ₊½(t, i, j, k, p) = ϕ_kedge(ζ̂z_val_and_ζderivs, t, i, j, k, p)

  return (;
    ξ̂x,
    η̂x,
    ζ̂x,
    ξ̂y,
    η̂y,
    ζ̂y,
    ξ̂z,
    η̂z,
    ζ̂z,
    ξ̂xᵢ₊½,
    η̂xᵢ₊½,
    ζ̂xᵢ₊½,
    ξ̂yᵢ₊½,
    η̂yᵢ₊½,
    ζ̂yᵢ₊½,
    ξ̂zᵢ₊½,
    η̂zᵢ₊½,
    ζ̂zᵢ₊½,
    ξ̂xⱼ₊½,
    η̂xⱼ₊½,
    ζ̂xⱼ₊½,
    ξ̂yⱼ₊½,
    η̂yⱼ₊½,
    ζ̂yⱼ₊½,
    ξ̂zⱼ₊½,
    η̂zⱼ₊½,
    ζ̂zⱼ₊½,
    ξ̂xₖ₊½,
    η̂xₖ₊½,
    ζ̂xₖ₊½,
    ξ̂yₖ₊½,
    η̂yₖ₊½,
    ζ̂yₖ₊½,
    ξ̂zₖ₊½,
    η̂zₖ₊½,
    ζ̂zₖ₊½,
  )
end

xmap, ymap, zmap = wavy_mapping()
celldims = (40, 40, 1000)
params = wavy_params(celldims)
backend = AutoForwardDiff()

eterms = get_edge_terms(xmap, ymap, zmap, backend)
nothing

#   t, i, j, k = (0.0, 10, 11.5, 12.5)

#   # @code_warntype ξ̂x(t, i, j, k, params)
#   # @benchmark ξ̂x($t, $i, $j, $k, $params)
#   idx = (i, j, k)
#   axis = 1
# end

# @code_warntype ζ̂xⱼ₊½(t, i, j, k, params)

# @code_warntype ζ̂x_val_and_ηderivs(t, k, j, k, params)

# ξ̂x_val_and_ξderivs = ξ_derivs(ξ̂x, backend)
# ξ̂x_iphalf = ϕ_edge(ξ̂x_val_and_derivs, t, idx, params, axis)
# @benchmark ϕ_edge($ξ̂x_val_and_derivs, $t, $idx, $params, $axis)

domain = CartesianIndices(celldims)

metrics = zeros(size(domain));

function compute_metrics(metrics, domain)
  # @profview begin
  for idx in domain
    i, j, k = idx.I
    # ξ̂x_half = ϕ_edge(ξ̂x_val_and_derivs, t, idx.I, params, axis)
    # ξ̂x_half = ϕ_edge(ξ̂x_val_and_derivs, t, i, j, k, params)
    metrics[idx] = eterms.η̂zₖ₊½(t, i, j, k, params)
  end
  return nothing
end

@profview compute_metrics(metrics, domain)
# end

# ϕall(t, i, j, k, params)
# @code_warntype ϕ_ᵢ₊½_3d(ϕall, t, i, j, k, params)

# @benchmark ϕ_ᵢ₊½_3d($ϕall, $t, $i, $j, $k, $params)
# @profview for i in 1:100_000
#   ϕ_ᵢ₊½_3d(ϕall, t, i, j, k, params)
# end

# @benchmark ϕ_ᵢ₊½_3d($ϕall, $t, $i, $j, $k, $params)
# @benchmark ϕ_edge($ϕall, $t, $idx, $params, $axis)

# ϕ_ᵢ₊½_3d(ϕall, t, i, j, k, params)

# @code_warntype ϕ_edge(ϕall, t, idx, params, axis)