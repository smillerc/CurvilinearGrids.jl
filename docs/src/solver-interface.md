# Solver Interface

Finite-volume solvers should interact with mapped and discrete grids through
conservative volumes and face geometry. That keeps solver code independent of
whether the grid came from an analytic map or imported coordinates.

## Face Flux Geometry

For a face selected by `loc` at cell index `idx`:

```julia
geom = face_flux_geometry(grid, idx, loc)

S = geom.metric_vector
A = geom.area
n = geom.normal
q = geom.coordinate
```

Here `S` is the outward-oriented conserved face metric vector in the grid's
physical basis. `A = norm(S)` and `n = S / A` when the area is positive.

This is the object a Riemann solver needs:

```julia
metric_scale = norm(S)
normal = S / metric_scale
face_flux = metric_scale * flux(flux_scheme, left, right, thermo, normal)
```

Pressure-split fluxes should also use the metric vector directly for pressure
force components:

```julia
pressure_flux_momentum = pressure * S
```

## Cell Update Pattern

A mapped finite-volume residual typically has the form:

```math
\frac{dU}{dt}
= -\frac{1}{V}
  \sum_\alpha
  \left(F_{\alpha,+} - F_{\alpha,-}\right),
```

where `V = cellvolume(grid, I)` and each face flux is assembled from
`face_flux_geometry`.

Code structure:

```julia
V = cellvolume(grid, I)

geom_hi = face_flux_geometry(grid, I.I, :ihi)
geom_lo = face_flux_geometry(grid, I.I, :ilo)

F_hi = riemann_flux(left_hi, right_hi, geom_hi.metric_vector)
F_lo = riemann_flux(left_lo, right_lo, geom_lo.metric_vector)

rhs = -(F_hi - F_lo) / V
```

## Why Not Read Raw Metrics Everywhere?

Raw face metric storage is useful for low-level algorithms and diagnostics:

```julia
Ghat = face_metrics(grid)[axis].conserved[I].jacobian_matrix
```

But solvers also need:

- side orientation,
- coordinate-system measure factors,
- area,
- normalized face normal,
- consistent face coordinates.

`face_flux_geometry` packages these into one API and is shared by `MappedGrid`
and `DiscreteGrid`.

## Orthogonal Grids

Orthogonal grids use analytical geometry APIs:

```julia
cellvolume(grid, I)
face_area(grid, I, axis)
outward_face_normal(grid, I, loc)
face_coordinate(grid, I, loc)
```

They intentionally do not expose mapped/discrete metric tensor caches.
