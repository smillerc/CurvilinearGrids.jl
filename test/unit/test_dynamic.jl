
using CurvilinearGrids, UnPack, Printf

function wavy_grid!(nodes, params)
  @unpack x, y = nodes

  τ = params.t

  Lx = Ly = 12
  Ax = Ay = 1.5
  nxy = nyx = 4.0
  ω = 1

  xmin = -Lx / 2
  ymin = -Ly / 2

  Δx0 = Lx / ni
  Δy0 = Ly / nj

  for idx in CartesianIndices(x)
    i, j = idx.I

    x[i, j] = xmin + Δx0 * (i + Ax * sinpi(2 * ω * τ) * sin((nxy * pi * Δy0 * j / Ly)))
    y[i, j] = ymin + Δy0 * (j + Ay * sinpi(2 * ω * τ) * sin((nyx * pi * Δx0 * i / Lx)))
  end

  return nothing
end

begin
  ni = nj = 50
  x = zeros(ni, nj)
  y = zeros(ni, nj)
  nodes = (; x, y)
  wavy_grid!(nodes, (; t=0))

  nhalo = 4
  mesh = CurvilinearGrid2D(x, y, nhalo; is_static=false)
  nothing

  Δt = 0.025
  t_all = 0:Δt:1
  for (cyc, tau) in enumerate(t_all)
    params = (; Δt=Δt, t=tau)
    CurvilinearGrids.update!(mesh, wavy_grid!, params)
    fn = "mesh" * @sprintf("%04i", cyc)

    save_vtk(mesh, fn)
  end
end
