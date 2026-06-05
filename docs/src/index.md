# CurvilinearGrids.jl

`CurvilinearGrids.jl` provides structured grid types, coordinate-system traits, and
metric caches for finite-volume and finite-difference solvers on curvilinear
meshes. The package is built around one central contract:

> A solver should be able to ask the grid for conservative cell volumes and
> solver-facing face geometry without knowing whether coordinates came from an
> analytic map, imported nodes, or an orthogonal coordinate family.

The current unified grid model has three main grid families:

| Grid type | Geometry source | Intended use |
|---|---|---|
| `MappedGrid` | analytic mapping functions | smooth generated geometry, manufactured tests, moving maps |
| `DiscreteGrid` | user-provided node coordinates | CAD/imported/generated meshes where only nodes are available |
| `OrthogonalGrid` | coordinate lines and analytic metric factors | Cartesian, cylindrical, spherical, and axisymmetric orthogonal grids |

## What To Read

- [Getting Started](getting-started.md): construct mapped and discrete grids.
- [Metric Schemes and GCL](metric-schemes.md): conservative metrics, GCL, legacy stencils, and AD Thomas-Lombard metrics.
- [Unified Grid API](unified-grid-api.md): common accessors for coordinates, metrics, and cache lifecycle.
- [Solver Interface](solver-interface.md): how finite-volume solvers consume face geometry.
- [Mesh Quality and Testing](mesh-quality.md): practical expectations for CAD/imported grids.
- [API Reference](api-reference.md): exported constructors, traits, and geometry functions.

## Design Position

`MappedGrid + AD` is the highest-fidelity path when an analytic mapping exists.
`DiscreteGrid + conservative metric reconstruction` is the production path for
imported coordinates. These are complementary, not interchangeable. The bridge
between them is verification: GCL residuals, face metric vectors, face areas,
face normals, and cell volumes must all be tested on the geometry that solvers
actually consume.

Existing design RFCs live in the package `docs/` directory and remain useful
historical context for the unified grid design.
