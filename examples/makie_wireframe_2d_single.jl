using CairoMakie
using CurvilinearGrids

# You can swap CairoMakie with GLMakie/WGLMakie for interactive rendering.

const NHALO = 5
const CELLDIMS = (24, 18)

xmap(t, ξ, η, p) = ξ + p.αx * sin(2π * η / p.Lη)
ymap(t, ξ, η, p) = η + p.αy * sin(2π * ξ / p.Lξ)

params = (; αx=0.20, αy=0.12, Lξ=CELLDIMS[1], Lη=CELLDIMS[2])

mapped = MappedGrid(
  xmap,
  ymap,
  params,
  CELLDIMS,
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

x_nodes, y_nodes = coords(mapped)
discrete = DiscreteGrid(
  collect(x_nodes),
  collect(y_nodes),
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

fig = Figure(; size=(1200, 520))
ax1 = Axis(fig[1, 1]; title="MappedGrid (2D)", aspect=DataAspect())
ax2 = Axis(fig[1, 2]; title="DiscreteGrid (2D)", aspect=DataAspect())

mesh!(ax1, mapped; color=(:steelblue, 0.28), transparency=true)
mesh!(ax2, discrete; color=(:seagreen, 0.28), transparency=true)
wireframe!(ax1, mapped; color=:black, linewidth=1.1, marker=:circle, markersize=3)
wireframe!(ax2, discrete; color=:black, linewidth=1.1, marker=:circle, markersize=3)

fig[0, :] = Label(
  fig, "CurvilinearGrids Makie mesh + wireframe extension (2D single-block)"; fontsize=18
)

display(fig)
save(joinpath(@__DIR__, "makie_wireframe_2d_single.png"), fig)
