using CairoMakie
using CurvilinearGrids

const NHALO = 5
const CELLDIMS = (10, 8, 7)
const SHIFT_X = CELLDIMS[1]

xleft(t, ξ, η, ζ, p) = ξ + p.α * sin(2π * η / p.Lη) * sin(2π * ζ / p.Lζ)
yleft(t, ξ, η, ζ, p) = η + p.β * sin(2π * ζ / p.Lζ)
zleft(t, ξ, η, ζ, p) = ζ + p.γ * sin(2π * η / p.Lη)

xright(t, ξ, η, ζ, p) = ξ + p.shift + p.α * sin(2π * η / p.Lη) * sin(2π * ζ / p.Lζ)
yright(t, ξ, η, ζ, p) = η + p.β * sin(2π * ζ / p.Lζ)
zright(t, ξ, η, ζ, p) = ζ + p.γ * sin(2π * η / p.Lη)

left_params = (; α=0.14, β=0.11, γ=0.10, Lη=CELLDIMS[2], Lζ=CELLDIMS[3], shift=0.0)
right_params = (; α=0.14, β=0.11, γ=0.10, Lη=CELLDIMS[2], Lζ=CELLDIMS[3], shift=SHIFT_X)

mapped_left = MappedGrid(
  xleft,
  yleft,
  zleft,
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
  zright,
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

xl, yl, zl = coords(mapped_left)
xr, yr, zr = coords(mapped_right)

discrete_left = DiscreteGrid(
  collect(xl),
  collect(yl),
  collect(zl),
  NHALO;
  compute_metrics=false,
  cache_mode=:off,
  coordinate_system=CurvilinearCS(),
  basis=CartesianBasis(),
)

discrete_right = DiscreteGrid(
  collect(xr),
  collect(yr),
  collect(zr),
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

fig = Figure(; size=(1280, 560))
ax1 = Axis3(fig[1, 1]; title="MappedGrid MultiBlock (3D)")
ax2 = Axis3(fig[1, 2]; title="DiscreteGrid MultiBlock (3D)")

mesh!(ax1, mapped_mb; transparency=true)
mesh!(ax2, discrete_mb; transparency=true)
wireframe!(ax1, mapped_mb; linewidth=0.9)
wireframe!(ax2, discrete_mb; linewidth=0.9)

fig[0, :] = Label(
  fig,
  "CurvilinearGrids Makie mesh + wireframe extension (3D multi-block boundaries)";
  fontsize=16,
)

display(fig)
save(joinpath(@__DIR__, "makie_wireframe_3d_multiblock.png"), fig)
