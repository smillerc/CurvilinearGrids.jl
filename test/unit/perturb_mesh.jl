
function gaussian(x, x0, fwhm, p)
  σ = fwhm / (2sqrt(2log(2)))
  σ² = σ * σ
  xc = (x - x0)^2
  return exp(-(((xc / (2σ²)))^p))
end

function perturb_coords!(mesh, x_interface, λ, k)
  for idx in mesh.iterators.node.domain
    i, j = idx.I

    x = mesh.node_coordinates.x[i, j]
    y = mesh.node_coordinates.y[i, j]
    amp = 0.25

    x_pert = (
      amp * # amplitude
      cos(k * y) * # mode
      gaussian(
        x,
        x_interface,
        0.5λ, # fwhm
        1.0, # gaussian exponent p
      ) # decay away from interface
    )

    mesh.node_coordinates.x[i, j] += x_pert
  end

  CurvilinearGrids.update!(mesh; force=true)
  return nothing
end

@testset "2D Mesh Perturbation" begin
  λ = 4.5 # wavelength
  k = 2pi / λ # wavenumber
  x_interface = 5.0

  x0, x1 = (0.0, 10.0)
  y0, y1 = (0.0, 0.5λ)
  nhalo = 4

  mesh = RectlinearGrid((x0, y0), (x1, y1), (501, 101), nhalo; is_static=true)

  perturb_coords!(mesh, x_interface, λ, k)

  # save_vtk(mesh)

  function gcl(mesh, domain)
    ϵ = 5e-14
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

      # @show I₁, I₂
      I₁_passes = abs(I₁) < eps()
      I₂_passes = abs(I₂) < eps()
      if !(I₁_passes && I₂_passes)
        break
      end
    end
    @test I₁_passes
    @test I₂_passes
  end

  gcl(mesh, mesh.iterators.cell.domain)
end