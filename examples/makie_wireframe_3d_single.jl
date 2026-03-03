using CairoMakie
using CurvilinearGrids

const NHALO = 5
const CELLDIMS = (14, 12, 10)

xmap(t, ξ, η, ζ, p) = ξ + p.αx * sin(2π * η / p.Lη)
ymap(t, ξ, η, ζ, p) = η + p.αy * sin(2π * ζ / p.Lζ)
zmap(t, ξ, η, ζ, p) = ζ + p.αz * sin(2π * ξ / p.Lξ)

params = (; αx=0.15, αy=0.12, αz=0.10, Lξ=CELLDIMS[1], Lη=CELLDIMS[2], Lζ=CELLDIMS[3])

mapped = MappedGrid(
  xmap,
  ymap,
  zmap,
  params,
  CELLDIMS,
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

x_nodes, y_nodes, z_nodes = coords(mapped)
discrete = DiscreteGrid(
  collect(x_nodes),
  collect(y_nodes),
  collect(z_nodes),
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

fig = Figure(; size=(1280, 560))
ax1 = Axis3(fig[1, 1]; title="MappedGrid (3D)")
ax2 = Axis3(fig[1, 2]; title="DiscreteGrid (3D)")

mesh!(ax1, mapped; color=(:steelblue, 0.32), transparency=true)
mesh!(ax2, discrete; color=(:seagreen, 0.32), transparency=true)
wireframe!(ax1, mapped; color=:black, linewidth=0.9)
wireframe!(ax2, discrete; color=:black, linewidth=0.9)

fig[0, :] = Label(
  fig,
  "CurvilinearGrids Makie mesh + wireframe extension (3D single-block boundaries)";
  fontsize=16,
)

display(fig)
save(joinpath(@__DIR__, "makie_wireframe_3d_single.png"), fig)
