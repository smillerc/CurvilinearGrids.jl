
using CurvilinearGrids
using UnPack
using Test

function gaussian(x, x0, fwhm, p)
  σ = fwhm / (2sqrt(2log(2)))
  σ² = σ * σ
  xc = (x - x0)^2
  return exp(-(((xc / (2σ²)))^p))
end

function rect_grid(nx, ny)
  λ = 4.5 # wavelength
  k = 2pi / λ # wavenumber

  x0, x1 = (0.0, 10.0)
  y0, y1 = (0.0, 0.5λ)

  x_interface = 5.0

  y1d(η) = @. y0 + (y1 - y0) * ((η - 1) / (ny - 1))
  x1d(ξ) = @. x0 + (x1 - x0) * ((ξ - 1) / (nx - 1))

  function x(ξ, η)
    xij = x1d(ξ)
    yij = y1d(η)

    amp = 0.25
    xpert = (
      amp * # amplitude
      cos(k * yij) * # mode
      gaussian(
        xij,
        x_interface,
        0.5λ, # fwhm
        1.0, # gaussian exponent p
      ) # decay away from interface
    )

    xij += xpert
    return xij
  end

  y(ξ, η) = y1d(η)

  return (x, y)
end

ni, nj = (101, 101)
nhalo = 4
x, y = rect_grid(ni, nj)
mesh = CurvilinearGrid2D(x, y, (ni, nj), nhalo)

domain = mesh.iterators.cell.domain

save_vtk(mesh)

# xn = mesh.node_coordinates.x[mesh.iterators.node.domain]
# yn = mesh.node_coordinates.y[mesh.iterators.node.domain]

# xc = mesh.centroid_coordinates.x[mesh.iterators.cell.domain]
# yc = mesh.centroid_coordinates.y[mesh.iterators.cell.domain]

# nodes = (xn, yn)
# centroids = (xc, yc)

# @unpack jacobian, ∂x_∂ξ, ∂x_∂η, ∂y_∂ξ, ∂y_∂η, ∂ξ_∂x, ∂η_∂x, ∂ξ_∂y, ∂η_∂y = CurvilinearGrids.FiniteDifferenceMetrics._compute_fd_metrics(
#   centroids, 4
# )
# # @unpack jacobian, ∂x_∂ξ, ∂x_∂η, ∂y_∂ξ, ∂y_∂η, ∂ξ_∂x, ∂η_∂x, ∂ξ_∂y, ∂η_∂y = CurvilinearGrids.FiniteDifferenceMetrics._compute_fd_metrics(
# #   nodes, 4
# # )

# ∂η_∂y
# ∂y_∂η
# begin
#   @test all(jacobian ≈ mesh.cell_center_metrics.J[domain])
#   @test all(∂ξ_∂x ≈ mesh.cell_center_metrics.ξ.x₁[domain])
#   @test all(∂ξ_∂y ≈ mesh.cell_center_metrics.ξ.x₂[domain])
#   @test all(∂η_∂x ≈ mesh.cell_center_metrics.η.x₁[domain])
#   @test all(∂η_∂y ≈ mesh.cell_center_metrics.η.x₂[domain])
# end
function gcl(mesh, domain)
  ϵ = 1e-14
  I₁_passes = true
  I₂_passes = true
  for idx in domain
    i, j = idx.I

    ξ̂_i₊½ = mesh.edge_metrics.i₊½.ξ̂[i, j]
    ξ̂_i₋½ = mesh.edge_metrics.i₊½.ξ̂[i - 1, j]
    η̂_j₊½ = mesh.edge_metrics.j₊½.η̂[i, j]
    η̂_j₋½ = mesh.edge_metrics.j₊½.η̂[i, j - 1]

    I₁ = (ξ̂_i₊½.x₁ - ξ̂_i₋½.x₁) + (η̂_j₊½.x₁ - η̂_j₋½.x₁)
    I₂ = (ξ̂_i₊½.x₂ - ξ̂_i₋½.x₂) + (η̂_j₊½.x₂ - η̂_j₋½.x₂)

    I₁ = I₁ * (abs(I₁) >= ϵ)
    I₂ = I₂ * (abs(I₂) >= ϵ)

    @show I₁, I₂
    I₁_passes = abs(I₁) < eps()
    I₂_passes = abs(I₂) < eps()
    if !(I₁_passes && I₂_passes)
      break
    end
  end
  @test I₁_passes
  @test I₂_passes
end

gcl(mesh, domain)
nothing