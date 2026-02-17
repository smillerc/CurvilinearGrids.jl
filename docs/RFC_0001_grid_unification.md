# RFC 0001: Unified Grid Type System and Traits

- Status: Draft
- Authors: CurvilinearGrids maintainers
- Depends on: `docs/RFC_0000_baseline.md`

## 1. Summary

This RFC proposes consolidating the grid model into exactly three grid types:

- `MappedGrid`: coordinates from continuous mapping functions
- `DiscreteGrid`: coordinates from user arrays, with an internal linear interpolated mapping
- `OrthogonalGrid`: orthogonal-coordinate geometry with no mapping metrics

The RFC also introduces traits to define:

- Coordinate system (all grid types)
- Basis vectors (only `MappedGrid` and `DiscreteGrid`)

For `MappedGrid` and `DiscreteGrid`, cell and face metrics are cached independently and invalidated by grid updates.

## 2. Motivation

Current grid taxonomy has many concrete types with overlapping behavior. This creates duplicated constructors, metric pipelines, and update logic. A 3-type model plus traits will:

- Reduce API surface complexity
- Unify operator dispatch rules
- Decouple coordinate-system semantics from concrete storage type
- Make new coordinate systems and bases extensible without type explosion

## 3. Goals and Non-Goals

Goals:
- Single conceptual model for curvilinear and orthogonal grids
- Trait-driven extensibility for coordinate systems and basis vectors
- Explicit metric cache lifecycle for mapped/discrete grids
- Backward-compatible migration path via adapters/deprecations

Non-goals:
- Rewriting all numerical schemes in this RFC
- Immediate removal of existing types in one release

## 4. Proposed Types

### 4.1 `MappedGrid`

`MappedGrid` stores mapping functions and evaluates coordinates/metrics from those mappings.

Indicative shape:

```julia
struct MappedGrid{N,CS,BT,Map,State,Cache,Iter,Backend}
  mapping::Map                # e.g. NamedTuple of x/y/z maps
  state::State                # mapping parameters, time, etc.
  metric_cache::Cache         # cell + face metric cache (optional/lazy)
  iterators::Iter             # node/cell + halo domains
  backend::Backend
end
```

### 4.2 `DiscreteGrid`

`DiscreteGrid` stores node coordinates from user arrays and derives a linear interpolated mapping used for metric evaluation/update.

Indicative shape:

```julia
struct DiscreteGrid{N,CS,BT,Coords,Interp,Cache,Iter,Backend}
  node_coordinates::Coords
  interpolant::Interp         # reconstruction from discrete nodes
  metric_cache::Cache         # independently managed cell + face metric caches
  iterators::Iter
  backend::Backend
end
```

### 4.3 `OrthogonalGrid`

`OrthogonalGrid` stores coordinate lines and analytical geometric terms (cell volumes, face areas, scale factors if needed). No mapping-metric tensors (`J`, `ξx`, etc.) are required.

Indicative shape:

```julia
struct OrthogonalGrid{N,CS,Coords,Geom,Iter,Backend}
  coordinates::Coords         # coordinate vectors like (x, y, z), or (r, z) or (r, θ, ϕ)
  geometry_cache::Geom        # volumes/areas/scale factors
  iterators::Iter
  backend::Backend
end
```

## 5. Trait System

## 5.1 Coordinate System Trait (all grid types)

```julia
abstract type CoordinateSystemTrait end
struct CartesianCS <: CoordinateSystemTrait end
struct CylindricalCS <: CoordinateSystemTrait end
struct SphericalCS <: CoordinateSystemTrait end
struct AxisymmetricCS{Axis} <: CoordinateSystemTrait end
struct CurvilinearCS <: CoordinateSystemTrait end

coordinate_system(::Type{<:AnyGrid}) -> CoordinateSystemTrait
```

Purpose:
- Determines geometry semantics, IO labeling, and operator specializations.

## 5.2 Basis Vector Trait (mapped/discrete only)

```julia
abstract type BasisTrait end
struct CartesianBasis <: BasisTrait end
struct ContravariantBasis <: BasisTrait end
struct CovariantBasis <: BasisTrait end
struct SphericalBasis <: BasisTrait end

basis_trait(::Type{<:Union{MappedGrid,DiscreteGrid}}) -> BasisTrait
```

Purpose:
- Encodes which basis vectors metric operators assume.
- Controls metric tensor fields to compute/cache.

Constraint:
- `basis_trait` is undefined for `OrthogonalGrid`.

## 6. Metrics and Caching

### 6.1 Cache Model

`MappedGrid` and `DiscreteGrid` cache, independently:
- Cell metrics cache: Jacobian and cell-centered forward/inverse metric tensors
- Face metrics cache: per-face metrics (`i+1/2`, `j+1/2`, `k+1/2`) including normalized conservative hatted metrics (`ξ̂`, `η̂`, `ζ̂`) as in existing grid types

Cache modes:
- `:eager` compute at construction/update
- `:lazy` compute on first access
- `:off` do not store; compute on demand (debug/reference path)

### 6.2 Cache Invalidation Rules

- Any coordinate/mapping/state mutation invalidates both caches.
- Cell and face caches may be recomputed independently.
- For `DiscreteGrid`, interpolation is fixed to linear; changing interpolation scheme is out of scope.
- Read APIs must ensure cache validity (recompute if stale and enabled).

### 6.3 Orthogonal Exception

`OrthogonalGrid` does not expose mapping metrics. It may cache:
- Cell volumes
- Face areas
- Optional scale factors

## 7. Unified API (Target)

Construction:
- `MappedGrid(mapping, dims; kwargs...)`
- `DiscreteGrid(node_coords; kwargs...)`
- `OrthogonalGrid(coords; coordinate_system=..., kwargs...)`

Common geometry:
- `coords(grid)`, `coord(grid, idx)`
- `centroids(grid)`, `centroid(grid, idx)`
- `cellvolume(grid, idx)`, `cellvolumes(grid)`

Mapped/discrete only:
- `cell_metrics(grid; refresh=false)`
- `face_metrics(grid; refresh=false)`
- `jacobian_matrix(grid, idx)`

Lifecycle:
- `update!(grid, args...)`
- `invalidate_cell_metrics!(grid)`
- `invalidate_face_metrics!(grid)`
- `refresh_cell_metrics!(grid; include_halo_region=false)`
- `refresh_face_metrics!(grid; include_halo_region=false)`

Trait queries:
- `coordinate_system(typeof(grid))`
- `basis_trait(typeof(grid))` (mapped/discrete)

## 8. Dispatch Strategy

Use trait-based internal dispatch rather than many concrete type branches:

```julia
_divergence(grid, field, idx, ::CartesianCS)
_divergence(grid, field, idx, ::SphericalCS)

_metric_kernel(grid, idx, ::ContravariantBasis)
_metric_kernel(grid, idx, ::SphericalBasis)
```

Primary method signatures stay grid-centric; traits are resolved internally.

## 9. Migration Plan

Phase 1 (additive):
- Introduce new types and traits behind stable constructors.
- Provide adapters from legacy grids to `MappedGrid`/`DiscreteGrid`/`OrthogonalGrid`.

Phase 2 (dual-path):
- Operators and IO accept both legacy and new grids.
- Emit deprecation warnings for legacy constructor entry points.

Phase 3 (cleanup):
- Remove legacy concrete types after one major release cycle.
- Keep compatibility shims for serialized formats where feasible.

## 10. Backward Compatibility

- Legacy constructors (`CurvilinearGrid*`, `RectilinearGrid*`, `UniformGrid*`, etc.) become wrappers returning new types.
- Existing top-level APIs (`coord`, `centroid`, `save_vtk`, `write_coordinates`) continue working through generic interfaces.
- Serialized `grid_type` strings in HDF5 should include version metadata to disambiguate old/new layouts.

## 11. Open Design Questions

- Should cache store both forward and inverse metrics by default, or compute inverse lazily?
- Do basis traits belong to type parameters or runtime trait functions only?
- How much orthogonal geometry (areas/volumes/scale factors) should be required vs lazily derived?
- Should `MappedGrid` state include time explicitly or be fully user-managed?

## 12. Acceptance Criteria

This RFC is accepted when:
- Three grid types are implemented and documented.
- Coordinate-system trait exists and is used across operators and IO.
- Basis trait exists and is used in mapped/discrete metric pathways.
- Metric cache lifecycle (`invalidate`/`refresh`) is implemented and tested.
- Legacy constructor wrappers pass existing core tests.

## 13. Minimal Test Matrix (for implementation RFC follow-up)

- Construction:
  - 1D/2D/3D for each new grid type where applicable
- Traits:
  - Coordinate trait correctness per constructor path
  - Basis trait correctness per mapped/discrete configuration
- Metrics:
  - Cache stale/refresh correctness
  - Consistency of cell vs face metric updates
- Operators:
  - Representative gradient/divergence on at least Cartesian and spherical systems
- Compatibility:
  - Legacy constructor wrappers produce equivalent geometry and metrics within tolerance
