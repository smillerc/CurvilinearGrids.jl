
using CurvilinearGrids, UnPack, Printf, Glob

function wavy_grid!(nodes, params)
  @unpack x, y = nodes
  @unpack t, Ax, Ay, Az, Lx, Ly, Lz, nxy, nyx, nzx, nzy, nxz, nyz, ω = params

  xmin = -Lx / 2
  ymin = -Ly / 2
  zmin = -Lz / 2

  Δx0 = Lx / ni
  Δy0 = Ly / nj
  Δz0 = Lz / nk

  for idx in CartesianIndices(x)
    i, j, k = idx.I

    #! format: off
    x[idx] = xmin + Δx0 * (i + Ax * sinpi(2 * ω * t) * sinpi((nxy * Δy0 * j / Ly)) * sinpi((nxz * Δz0 * k / Lz)))
    y[idx] = ymin + Δy0 * (j + Ay * sinpi(2 * ω * t) * sinpi((nyx * Δx0 * i / Lx)) * sinpi((nyz * Δz0 * k / Lz)))
    z[idx] = zmin + Δz0 * (k + Az * sinpi(2 * ω * t) * sinpi((nzx * Δx0 * i / Lx)) * sinpi((nzy * Δy0 * j / Ly)))
    #! format: on

  end

  return nothing
end

function grid_vel!(velocity, params)
  @unpack x, y = velocity
  @unpack t, Lx, Ly, nxy, nyx, ω = params

  Δx0 = Lx / ni
  Δy0 = Ly / nj

  for idx in CartesianIndices(x)
    i, j, k = idx.I

    #! format: off
    x[idx] = 2pi * ω * Δx0 * Ax * cospi(2 * ω * t) * sinpi((nxy * Δy0 * j / Ly)) * sinpi((nxz * Δz0 * k / Lz))
    y[idx] = 2pi * ω * Δy0 * Ay * cospi(2 * ω * t) * sinpi((nyx * Δx0 * i / Lx)) * sinpi((nyz * Δz0 * k / Lz))
    z[idx] = 2pi * ω * Δz0 * Az * cospi(2 * ω * t) * sinpi((nzx * Δx0 * i / Lx)) * sinpi((nzy * Δy0 * j / Ly))
    #! format: on

  end

  return nothing
end

begin
  rm.(glob("*.vts"))

  ni = nj = 50
  x = zeros(ni, nj)
  y = zeros(ni, nj)
  nodes = (; x, y)
  nxy = nyx = nzx = nzy = nxz = nyz = 4.0
  grid_params = (; nxy, nyx, nzx, nzy, nxz, nyz, Lx=12, Ly=12, Lz=12, ω=1, Ax=1.5, Ay=1.5)
  params = (; Δt=0, t=0)
  wavy_grid!(nodes, merge(grid_params, params))

  nhalo = 4
  mesh = CurvilinearGrid2D(
    x, y, nhalo; coord_update_function=wavy_grid!, coord_velocity_function=grid_vel!
  )

  Δt = 0.05
  t_all = 0:Δt:1
  for (cyc, tau) in enumerate(t_all)
    t_params = (; Δt=Δt, t=tau)
    p = merge(t_params, grid_params)

    CurvilinearGrids.update!(mesh, p)
    fn = "mesh" * @sprintf("%04i", cyc)

    save_vtk(mesh, fn)
  end
end
