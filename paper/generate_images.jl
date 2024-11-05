using CairoMakie
using CurvilinearGrids

function wavy_grid_1(nx, ny, nhalo=1)
  x2d = zeros(nx, ny)
  y2d = zeros(nx, ny)

  x1d = range(0, 1; length=nx)
  y1d = range(0, 1; length=ny)
  a0 = 0.1
  for I in CartesianIndices(x2d)
    i, j = I.I

    x = x1d[i]
    y = y1d[j]

    x2d[i, j] = x + a0 * sinpi(2x) * sinpi(2y)
    y2d[i, j] = y + a0 * sinpi(2x) * sinpi(2y)

    # x2d[i, j] = ((i - 1) / (nx - 1)) + a0 * sinpi(2x) * sinpi(2y)
    # y2d[i, j] = ((j - 1) / (ny - 1)) + a0 * sinpi(2x) * sinpi(2y)
  end

  mesh = CurvilinearGrids.CurvilinearGrid2D(x2d, y2d, nhalo)
  return mesh
end

function wavy_grid_2(ni, nj, nhalo=1)
  Lx = 12
  Ly = 12
  n_xy = 6
  n_yx = 6

  xmin = -Lx / 2
  ymin = -Ly / 2

  Δx0 = Lx / (ni - 1)
  Δy0 = Ly / (nj - 1)

  Ax = 0.4 / Δx0
  Ay = 0.8 / Δy0

  x = zeros(ni, nj)
  y = zeros(ni, nj)
  for j in 1:nj
    for i in 1:ni
      x[i, j] = xmin + Δx0 * ((i - 1) + Ax * sinpi((n_xy * (j - 1) * Δy0) / Ly))
      y[i, j] = ymin + Δy0 * ((j - 1) + Ay * sinpi((n_yx * (i - 1) * Δx0) / Lx))
    end
  end

  return CurvilinearGrid2D(x, y, nhalo)
end

function uniform_grid(nx, ny, nhalo=1)
  x0, x1 = (0, 1)
  y0, y1 = (0, 1)

  mesh = CurvilinearGrids.RectlinearGrid((x0, y0), (x1, y1), (nx, ny), nhalo)
  return mesh
end


nx = 30
ny = 40
grid = wavy_grid_2(nx, ny)

f = Figure(; size=(800, 500))
ax1 = Axis(
  f[1, 1]; aspect=1, xlabel="x", ylabel="y", xgridvisible=false, ygridvisible=false
)
ax2 = Axis(
  f[1, 2]; aspect=1, xlabel="ξ", ylabel="η", xgridvisible=false, ygridvisible=false
)

x, y = coords(grid)
lw = 0.5
ms = 3
for (_x, _y) in zip(eachcol(x), eachcol(y))
  scatterlines!(ax1, _x, _y; color=:black, linewidth=lw, markersize=ms)
end

for (_x, _y) in zip(eachrow(x), eachrow(y))
  scatterlines!(ax1, _x, _y; color=:black, linewidth=lw, markersize=ms)
end

# uniform grid
for i in 1:nx
  scatterlines!(ax2, i * ones(ny), 1:ny; color=:black, linewidth=lw, markersize=ms)
end

for i in 1:ny
  scatterlines!(ax2, 1:nx, i * ones(nx); color=:black, linewidth=lw, markersize=ms)
end

display(f)
save("mesh.png", f)

