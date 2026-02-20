
function gaussian(x, x0, fwhm, p)
  σ = fwhm / (2sqrt(2log(2)))
  σ² = σ * σ
  xc = (x - x0)^2
  return exp(-(((xc / (2σ²)))^p))
end

function perturb_coords!(mesh, x_interface, λ, k)
  xcoords = Array(mesh.node_coordinates.x)
  ycoords = Array(mesh.node_coordinates.y)

  for idx in mesh.iterators.node.domain
    i, j = idx.I

    x = xcoords[i, j]
    y = ycoords[i, j]
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

    xcoords[i, j] += x_pert
  end

  return DiscreteGrid(
    xcoords, ycoords, mesh.discretization_scheme_name; halo_coords_included=true
  )
end

@testset "2D Mesh Perturbation" begin
  λ = 4.5 # wavelength
  k = 2pi / λ # wavenumber
  x_interface = 5.0

  x0, x1 = (0.0, 10.0)
  y0, y1 = (0.0, 0.5λ)

  xnodes = collect(range(x0, x1, 502))
  ynodes = collect(range(y0, y1, 102))
  x = [xnodes[i] for i in eachindex(xnodes), j in eachindex(ynodes)]
  y = [ynodes[j] for i in eachindex(xnodes), j in eachindex(ynodes)]
  mesh = DiscreteGrid(x, y, :meg6_symmetric)

  mesh = perturb_coords!(mesh, x_interface, λ, k)

  # save_vtk(mesh)

  gcl_identities, max_vals = gcl(face_metrics(mesh), mesh.iterators.cell.domain, eps())
  @test_broken all(gcl_identities)
end
