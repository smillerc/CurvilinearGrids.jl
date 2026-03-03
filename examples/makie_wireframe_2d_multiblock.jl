using CairoMakie
using CurvilinearGrids

# Two abutting blocks in i-direction: (left :ihi) <-> (right :ilo)

const NHALO = 5
const CELLDIMS = (20, 16)
const SHIFT_X = CELLDIMS[1]

xleft(t, ξ, η, p) = ξ + p.α * sin(2π * η / p.Lη)
yleft(t, ξ, η, p) = η

xright(t, ξ, η, p) = ξ + p.shift + p.α * sin(2π * η / p.Lη)
yright(t, ξ, η, p) = η

left_params = (; α=0.18, Lη=CELLDIMS[2], shift=0.0)
right_params = (; α=0.18, Lη=CELLDIMS[2], shift=SHIFT_X)

mapped_left = MappedGrid(
  xleft,
  yleft,
  left_params,
  CELLDIMS,
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

mapped_right = MappedGrid(
  xright,
  yright,
  right_params,
  CELLDIMS,
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

mapped_interfaces = ((mapped_left, :ihi) => (mapped_right, :ilo),)
mapped_mb = MultiBlockMesh((mapped_left, mapped_right), mapped_interfaces; tolerance=1e-10)

xl, yl = coords(mapped_left)
xr, yr = coords(mapped_right)

discrete_left = DiscreteGrid(
  collect(xl),
  collect(yl),
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

discrete_right = DiscreteGrid(
  collect(xr),
  collect(yr),
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

discrete_interfaces = ((discrete_left, :ihi) => (discrete_right, :ilo),)
discrete_mb = MultiBlockMesh(
  (discrete_left, discrete_right), discrete_interfaces; tolerance=1e-10
)

fig = Figure(; size=(1200, 520))
ax1 = Axis(fig[1, 1]; title="MappedGrid MultiBlock (2D)", aspect=DataAspect())
ax2 = Axis(fig[1, 2]; title="DiscreteGrid MultiBlock (2D)", aspect=DataAspect())

mesh!(ax1, mapped_mb; transparency=true)
mesh!(ax2, discrete_mb; transparency=true)
wireframe!(ax1, mapped_mb; linewidth=1.2, marker=:circle, markersize=2)
wireframe!(ax2, discrete_mb; linewidth=1.2, marker=:circle, markersize=2)

fig[0, :] = Label(
  fig, "CurvilinearGrids Makie mesh + wireframe extension (2D multi-block)"; fontsize=18
)

display(fig)
save(joinpath(@__DIR__, "makie_wireframe_2d_multiblock.png"), fig)
