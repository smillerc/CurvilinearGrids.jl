# RFC: CurvilinearGrids Specification Sheet

- Status: Draft
- Package: `CurvilinearGrids.jl`
- Version Baseline: `0.11.2`

## 1. Purpose

`CurvilinearGrids.jl` provides structured grid types and metric machinery for PDE discretizations on curvilinear meshes, with emphasis on conservative metric evaluation and halo-aware indexing.

## 2. Scope

This RFC covers:
- Discrete grid types (`CurvilinearGrid*`, `RectilinearGrid*`, `UniformGrid*`, axisymmetric/cylindrical/spherical variants)
- Continuous mapping-based grids (`ContinuousCurvilinearGrid*`)
- Metric update pipeline (forward + conservative inverse metrics + edge interpolation)
- IO (`save_vtk`, `write_coordinates`, `read_coordinates`)
- Resolution remapping helpers
- Spherical finite-volume differential operators (`gradient`, `divergence`, `derivative`)

Out of scope:
- Solver implementation
- Full boundary-condition framework
- Curl operators (present but disabled)

## 3. Design Goals

- Maintain conservative inverse metrics suitable for geometric conservation law (GCL) sensitive schemes.
- Separate node, cell-center, and edge representations.
- Keep domain iteration halo-aware.
- Support CPU/GPU-compatible arrays via `KernelAbstractions`.
- Support specialized orthogonal coordinate families where geometric terms are analytical.

## 4. Core Concepts

- Node coordinates: vertex locations.
- Centroid coordinates: per-cell centers.
- Cell metrics: Jacobian and forward/inverse metric tensors at cell centers.
- Edge metrics: interpolated metric tensors on `i+1/2`, `j+1/2`, `k+1/2`.
- Halo region: ghost-padding around the physical domain used for stencils and coupling.

## 5. Grid Families

### 5.1 Mapped/discrete grids

- 1D: `CurvilinearGrid1D`, `UniformGrid1D`, `CylindricalGrid1D`, `SphericalGrid1D`
- 2D: `CurvilinearGrid2D`, `RectilinearGrid2D`, `UniformGrid2D`, `AxisymmetricGrid2D`
- 3D: `CurvilinearGrid3D`, `RectilinearGrid3D`, `UniformGrid3D`, `SphericalBasisCurvilinearGrid3D`

Behavior:
- Construct from coordinate arrays/vectors and a discretization scheme symbol.
- Compute centroid coordinates.
- Compute cell and edge metrics (unless `empty_metrics=true`).
- Allow `update!(mesh, force, include_halo_region)` for dynamic meshes.

### 5.2 Orthogonal analytic grids

- `CartesianOrthogonalGrid1D`
- `CylindricalOrthogonalGrid1D`
- `SphericalOrthogonalGrid1D`
- `AxisymmetricOrthogonalGrid2D`
- `SphericalGrid3D`

Behavior:
- Geometric volumes/areas computed analytically from 1D coordinates.
- Used with dedicated spherical operators in `src/operators`.
- `SphericalGrid3D` stores native spherical node/centroid coordinates plus Cartesian node coordinates.

### 5.3 Continuous mapping grids

- `ContinuousCurvilinearGrid1D/2D/3D`

Behavior:
- Construct from mapping functions (e.g., `x(t, ξ, η, params)`).
- Auto-differentiate metric terms via `DifferentiationInterface` backend.
- `update!(mesh, t, params)` recomputes coordinates + metrics.

## 6. Discretization and Metric Schemes

- Main scheme family: `MonotoneExplicitGradientScheme` (MEG).
- Declared derivative orders: 2/4/6 (`:MEG2`, `:MEG4`, `:MEG6`, plus `:MEG6_SYMMETRIC`).
- Halo requirements (from scheme lookup):
  - `MEG2 -> nhalo=2`
  - `MEG4 -> nhalo=3`
  - `MEG6/MEG6_SYMMETRIC -> nhalo=5`
- Metric update sequence:
  1. forward derivatives (`x_ξ`, etc.)
  2. conservative inverse metrics (`ξ_x`, etc.), optionally symmetric conservative form
  3. edge interpolation of metrics

## 7. Public API Surface (Current)

Grid lifecycle and geometry:
- `update!`
- `coord`, `coords`
- `centroid`, `centroids`, `cartesian_centroid`
- `cellsize`, `cellsize_withhalo`
- `cellvolume`, `cellvolumes`
- `radius`, `centroid_radius`, `centroid_radii`
- `jacobian_matrix`
- `forward_cell_metrics`, `inverse_cell_metrics`

Constructors/helpers:
- `rectilinear_grid`
- `rectilinear_cylindrical_grid`
- `rectilinear_spherical_grid`
- `axisymmetric_rectilinear_grid`
- `rtheta_grid`, `axisymmetric_rtheta_grid`
- `rthetaphi_grid`

Operators (currently spherical-focused):
- `cell_center_derivative`, `edge_derivative`
- `cell_center_divergence`, `edge_divergence`
- `cell_center_gradient`, `edge_gradient`

IO:
- `save_vtk`
- `write_coordinates`, `read_coordinates`

Remapping:
- `change_resolution`
- `scale_resolution`
- `remap_cell_data`

## 8. Invariants and Contracts

- Coordinate vectors for rectilinear/uniform constructors must be strictly increasing.
- Jacobian must remain finite and positive in valid meshes; update paths perform metric validity checks for most mapped grids.
- Halo-aware index semantics:
  - `coord/centroid` accept halo-aware indices.
  - `coords/centroids` return interior-domain views (no halo).
- Static meshes:
  - `is_static=true` prevents accidental recompute unless `force=true`.
- Edge metrics are defined on expanded edge domains and are expected to satisfy GCL consistency checks.

## 9. Data Layout

- Metrics and coordinates are stored predominantly in `StructArray`-style field-of-arrays layout.
- Cell metric components include:
  - `J`
  - forward basis (`x₁`, `x₂`, `x₃` components wrt `ξ/η/ζ`)
  - inverse basis (`ξ`, `η`, `ζ` wrt physical axes)
  - normalized inverse basis (`ξ̂`, `η̂`, `ζ̂`)
- Edge metrics keyed by edge family (`i₊½`, `j₊½`, `k₊½`).

## 10. Testing Coverage Snapshot

The test suite includes:
- 1D/2D/3D mapped grids
- Continuous grid variants
- Axisymmetric and wall-distance utilities
- Discretization correctness checks for MEG derivatives
- HDF5 read/write round-trips for several grid families
- Rectilinear array behavior and perturbation scenarios

Notably disabled or partial:
- Orthogonal reduced grids tests are commented in `test/runtests.jl`.
- Curl operator tests absent (curl implementation disabled).

## 11. Known Gaps / Risks (Observed in Current Code)

- `src/io/to_h5.jl`: `read_SphericalGrid3D` returns `collect(x)` where `x` is undefined; expected `collect(r)`.
- `src/grids/curvilinear_mapped_grids/3d.jl`: `jacobian_matrix(mesh::CurvilinearGrid3D, ...)` appears to reuse `x₂`/`η` components for several matrix entries instead of using full `(x,y,z)` metric components.
- `src/grids/curvilinear_mapped_grids/3d_spherical_basis.jl`: analogous `jacobian_matrix` component reuse issue.
- `src/remap.jl`: `change_resolution(mesh::CurvilinearGrid1D, ni)` references `nj` and treats 1D interpolation as 2D.
- `src/remap.jl`: `scale_resolution(mesh::CurvilinearGrid1D, α)` references undefined `β`, `new_y`, and returns a 2D grid.
- Error message in scheme lookup states only `:MEG6` supported, while code supports `:MEG2`, `:MEG4`, `:MEG6`, and `:MEG6_SYMMETRIC`.

These should be considered implementation defects relative to the intended API contract.

## 12. Recommended Next RFCs

- RFC-A: Formal scheme support matrix and deprecation policy (`MEG2/4/6` status).
- RFC-B: Unified operator interface across mapped and orthogonal families.
- RFC-C: IO format versioning and compatibility guarantees.
- RFC-D: Validation harness for GCL, Jacobian positivity, and halo contract compliance.

