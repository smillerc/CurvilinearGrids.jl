# Changelog

All notable changes to this project are documented in this file.

## [Unreleased]

### Added
- Unified grid API centered on `MappedGrid`, `DiscreteGrid`, and `OrthogonalGrid`.
- Trait-based coordinate/basis dispatch (`coordinate_system`, `basis_trait`) for unified grids.
- Public solver-facing basis/geometry helpers:
  - `basis_transfer_matrix`
  - `face_coordinate`
  - `FaceFluxGeometry`
  - `face_flux_geometry`
- Continuous-coordinate evaluation for mapped/discrete grids in:
  - `coord`
  - `jacobian_matrix`
  - `forward_cell_metrics`
  - `inverse_cell_metrics`
  - `cellvolume`
- Inverse mapping API from physical to computational space:
  - `computational_coordinate`
  - `InverseCoordinateResult`
- Optional metric cache controls for mapped/discrete grids with `cache_mode = :eager | :lazy | :off`.
- Full metric-storage disable path via `compute_metrics=false, cache_mode=:off`.
- Independent cache invalidation/refresh API:
  - `invalidate_cell_metrics!`, `invalidate_face_metrics!`
  - `refresh_cell_metrics!`, `refresh_face_metrics!`
- Multi-block interface system:
  - `MultiBlockMesh`
  - `BlockFace`
  - `BlockInterface`
  - `validate_multiblock!`
  - `build_interface_caches!`, `invalidate_interface_caches!`
  - `exchange_interface!`, `exchange_all_interfaces!`
  - `interface_index_plan`
  - multi-block `computational_coordinate` transfer helpers
- Public scalar index-planning API for oriented multi-block interfaces:
  `interface_index_plan` returns ordered source-to-target Cartesian index pairs
  for both interface directions, including permutation, flip, halo-depth, and
  optional valid-domain filtering.
- New orthogonal-grid constructors for Cartesian grids:
  - `CartesianOrthogonalGrid2D`
  - `CartesianOrthogonalGrid3D`
- Public geometry/metric helpers for unified operators:
  - `cell_jacobian(grid, idx)`
  - `face_metric_coefficient(grid, dim, idx)`
  - axis-indexed `face_area(grid, dim, idx)` for orthogonal grids
- Orthogonal-grid `centroid`/`centroids` are part of the unified public
  solver-facing API, returning native physical coordinates as `SVector`s and
  coordinate-array views.

### Changed
- `DiscreteGrid` metric Jacobians are evaluated via `Interpolations.gradient` for interpolation-consistent derivatives.
- Unified-grid cell volume logic is now coordinate-system and basis aware.
- Multiblock basis exchange now goes through the same public
  `basis_transfer_matrix` API used by solver code.
- Multi-block scalar interface exchange now shares its oriented layer-mapping
  logic with `interface_index_plan`, so solver and MPI layers can reuse the same
  index plan instead of duplicating flip/permutation handling.
- HDF5 reconstruction for unified grids now supports `compute_metrics` and `cache_mode` options in `read_coordinates`.
- Unified-grid internals were refactored toward immutable, trait-driven storage and clearer cache/state behavior.
- `nhalo` handling was cleaned up in unified-grid construction paths.
- Orthogonal-grid geometry access now dispatches on generic `OrthogonalGrid{N}` for:
  - `coords`, `coord`
  - `centroids`, `centroid`
  - `cellvolume`
- Mapped/discrete face-flux geometry now has a public adapter over conserved
  inverse face metrics, with `face_flux_geometry(...).metric_vector` defined as
  the outward-oriented active conserved row in the grid's physical basis.
- Orthogonal-grid constructor internals now consistently honor `halo_coords_included` (including `nhalo == 0` face-domain paths).
- Backward-compatible identity constructor `OrthogonalGrid(grid::OrthogonalGrid)` is available for generic call sites.

### Fixed
- Missing exports for multi-block APIs from the top-level package.
- I/O compatibility paths for mapped/discrete grid serialization/deserialization.
- Spherical 3D orthogonal-grid centroid and cell-volume computations now use
  algebraically equivalent, numerically stabler radial and polar-angle
  formulas.

### Removed
- Legacy orthogonal-grid concrete type aliases were removed:
  - `CartesianOrthogonalGrid1D`
  - `CylindricalOrthogonalGrid1D`
  - `SphericalOrthogonalGrid1D`
  - `AxisymmetricOrthogonalGrid2D`
  - `SphericalGrid3D` (as a type alias)
- Constructor functions with these names remain available and return `OrthogonalGrid`.

### Refactored
- Multi-block functionality moved into a dedicated top-level `MultiBlockMeshes` module.
- Cubed-sphere helper moved to `mesh_functions`.
- Ongoing readability and test-coverage updates for unified-grid and multiblock workflows.
- Repository-wide JuliaFormatter normalization across unified-grid, orthogonal-grid, and operator sources.

### Documentation
- `README.md` was rewritten around the unified-grid API (`MappedGrid`/`DiscreteGrid`/`OrthogonalGrid`), cache controls, inverse mapping, and multi-block usage.
- Expanded the public documentation for solver-facing geometry and basis APIs,
  including the distinction between `outward_face_normal` and
  `face_flux_geometry(...).normal`.
- Added `docs/RFC_0002_multiblock_interfaces.md` describing the multi-block interface model and phased implementation plan.
- Expanded `docs/RFC_1111_coding_conventions.md` with explicit import rules and HPC-oriented coding-style guidance.

## Previous releases

- Releases through `v0.11.2` are available via git tags.
