
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
  xξ(t, i, j, k, p) = derivative(ξ -> x(t, ξ, j, k, p), backend, i)

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