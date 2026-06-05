# Unified Grid API

The unified API is organized around geometry source, coordinate semantics, and
solver-facing metric access.

## Grid Families

### `MappedGrid`

`MappedGrid` is for analytic coordinate maps:

```julia
MappedGrid(xmap, params, celldims, nhalo; kwargs...)
MappedGrid(xmap, ymap, params, celldims, nhalo; kwargs...)
MappedGrid(xmap, ymap, zmap, params, celldims, nhalo; kwargs...)
```

Use it when mapping functions are part of the model and can be evaluated at
arbitrary computational locations.

### `DiscreteGrid`

`DiscreteGrid` is for user-provided node arrays:

```julia
DiscreteGrid(xnodes, nhalo; kwargs...)
DiscreteGrid(xnodes, ynodes, nhalo; kwargs...)
DiscreteGrid(xnodes, ynodes, znodes, nhalo; kwargs...)
```

Use it for imported meshes, CAD-generated meshes, or any coordinate set without
a trusted analytic map.

### `OrthogonalGrid`

`OrthogonalGrid` is for coordinate-line grids whose geometry has analytic
volumes and areas. Orthogonal grids do not expose mapping metric tensors through
`cell_metrics` and `face_metrics`.

## Common Traits

```julia
coordinate_system(grid)
basis_trait(grid)
```

Supported coordinate-system traits include:

```julia
CurvilinearCS()
CartesianCS()
CylindricalCS()
SphericalCS()
AxisymmetricCS{:x}()
AxisymmetricCS{:y}()
```

Mapped and discrete grids use basis traits:

```julia
CartesianBasis()
SphericalBasis()
```

## Coordinates

```julia
coords(grid)
coord(grid, idx)
cartesian_coordinates(grid)

centroids(grid)
centroid(grid, idx)
cartesian_centroids(grid)
cartesian_centroid(grid, idx)

face_coordinates(grid)
face_coordinate(grid, idx, loc)
```

Face locations use symbols such as `:ilo`, `:ihi`, `:jlo`, `:jhi`, `:klo`, and
`:khi`.

## Metrics and Volumes

Mapped and discrete grids expose cell and face metric storage:

```julia
cell_metrics(grid)
face_metrics(grid)

forward_cell_metrics(grid, idx)
inverse_cell_metrics(grid, idx)

jacobian_matrix(grid, idx)
cellvolume(grid, idx)
```

`face_metrics(grid)[axis].conserved[I]` stores conserved metric matrices. A
solver usually should not consume this raw storage directly unless it is
implementing a low-level metric operation. Prefer `face_flux_geometry`.

## Face Geometry

```julia
geom = face_flux_geometry(grid, idx, :ihi)

geom.coordinate
geom.metric_vector
geom.area
geom.normal
```

`metric_vector` is the outward-oriented physical conserved face metric vector.
It includes coordinate-system measure factors, such as axisymmetric `2*pi*r` or
spherical basis factors.

## Cache Lifecycle

Mapped and discrete grids can cache metrics eagerly, lazily, or not at all:

```julia
cache_mode=:eager
cache_mode=:lazy
cache_mode=:off
```

The public lifecycle hooks are:

```julia
invalidate_cell_metrics!(grid)
invalidate_face_metrics!(grid)
refresh_cell_metrics!(grid)
refresh_face_metrics!(grid)
update!(grid, args...)
```

Use explicit refresh/invalidation when coordinates or mapping parameters change
outside the normal constructor/update path.
