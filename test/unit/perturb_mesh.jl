
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

  mesh = CurvilinearGrid2D((x0, y0), (x1, y1), (501, 101), :meg6_symmetric; is_static=true)

  perturb_coords!(mesh, x_interface, λ, k)

  # save_vtk(mesh)

  gcl(mesh)
end
