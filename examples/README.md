# Makie Wireframe Examples

These scripts demonstrate the `CurvilinearGridsMakieExt` plotting recipes (`mesh` + `wireframe`) for:

- 2D single-block mapped and discrete grids
- 2D multi-block mapped and discrete grids
- 3D single-block mapped and discrete grids
- 3D multi-block mapped and discrete grids

All examples disable metric computation/storage for speed:

- `compute_metrics=false`
- `cache_mode=:off`

## Run

From the package root:

```bash
julia --project examples/makie_wireframe_2d_single.jl
julia --project examples/makie_wireframe_2d_multiblock.jl
julia --project examples/makie_wireframe_3d_single.jl
julia --project examples/makie_wireframe_3d_multiblock.jl
```

Each script overlays:

- `mesh!` for filled surfaces
- `wireframe!` for mesh lines

and saves a PNG next to the script.

## Backend choice

Scripts currently use `CairoMakie` for consistent file output. You can switch to
`GLMakie`/`WGLMakie` by replacing `using CairoMakie` in each script.
