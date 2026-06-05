# Metric Schemes and GCL

Curvilinear mapped solvers live or die by their metric scheme. The package
therefore distinguishes between pointwise coordinate derivatives, conservative
metric terms, and the solver-facing face geometry.

## Metric Notation

Let physical coordinates be generated from computational coordinates:

```math
x = x(\xi, \eta, \zeta).
```

The forward mapping Jacobian is

```math
F_{i\alpha} = \frac{\partial x_i}{\partial \xi_\alpha},
```

with determinant `J = det(F)`. The inverse metric tensor is

```math
G_{\alpha i} = \frac{\partial \xi_\alpha}{\partial x_i}.
```

Finite-volume solvers usually need the conserved inverse metric

```math
\hat{G}_{\alpha i} = J G_{\alpha i}.
```

On a face normal to computational axis `alpha`, the active row
`\hat{G}_{\alpha,:}` is the conservative face metric vector before
coordinate-system measure factors are applied.

## Geometric Conservation Law

For a stationary grid, the discrete GCL requires the discrete divergence of each
contravariant basis vector to vanish:

```math
\sum_\alpha \delta_\alpha \hat{G}_{\alpha i} = 0.
```

This is an algebraic statement about the metric discretization. It is necessary
for freestream preservation, but it is not sufficient to prove that two metric
schemes produce equivalent face geometry. A stronger check compares the actual
face metric vectors, face areas, normals, and coordinates consumed by a solver.

## Legacy Conservative Metric Stencils

The legacy discrete path computes conservative metrics from coordinate data using
discrete derivative and reconstruction operators. In 3D it follows the
Thomas-Lombard idea: write metric cofactors as differences of metric potentials
so that mixed derivatives cancel under compatible discrete operators.

One representative component is

```math
\hat{\xi}_x =
\partial_\zeta(y_\eta z) -
\partial_\eta(y_\zeta z).
```

The important property is not that this exact expression is magical; it is that
the same discrete derivative algebra appears in paired terms, so GCL cancellation
is preserved.

## Available Metric Scheme Selectors

Metric scheme selection is controlled by the `conserved_metric_scheme` keyword on
`MappedGrid` and `DiscreteGrid`.

```julia
grid = MappedGrid(
  xmap,
  ymap,
  params,
  (nx, ny),
  nhalo;
  conserved_metric_scheme=CurvatureCorrectedReconstruction(),
)
```

The available selectors are:

| Scheme | What it changes | Best use |
|---|---|---|
| `EndpointAverageReconstruction()` | Reconstructs face quantities from endpoint averages. | Debugging, low-complexity comparisons, very smooth simple grids. |
| `GradientCorrectedReconstruction()` | Uses endpoint values and first-derivative corrections. | Intermediate accuracy checks and sensitivity studies. |
| `CurvatureCorrectedReconstruction()` | Uses endpoint values plus first- and second-derivative corrections. | Default recommendation for mapped and discrete grids unless there is a specific reason to choose otherwise. |
| `ADThomasLombardMetric()` | In 3D on CPU, uses AD-evaluated Thomas-Lombard metric potentials to assemble conserved face metric rows. In 1D/2D it falls back to the direct metric form with a warning. | Analytic 3D `MappedGrid` cases where the map is known and high-quality face metrics are needed. |

The `diff_backend` keyword controls how mapping derivatives are evaluated inside
these paths. The default is `AutoForwardDiff()`.

```julia
grid = MappedGrid(
  xmap,
  ymap,
  zmap,
  params,
  (nx, ny, nz),
  nhalo;
  diff_backend=AutoForwardDiff(),
  conserved_metric_scheme=ADThomasLombardMetric(),
)
```

## Metric Storage Paths

Mapped and discrete grids store cell metrics and face metrics separately:

```julia
cell_metrics(grid)
face_metrics(grid)
```

This distinction matters. A face metric scheme does not necessarily use stored
cell metric arrays as its input.

| Storage | Contents | How it is used |
|---|---|---|
| `cell_metrics.forward` | cell-centered forward Jacobian `F` and determinant `J` | cell volumes, diagnostics, operators needing cell-centered geometry |
| `cell_metrics.inverse` | cell-centered inverse metric `G` and determinant `J` | cell-centered inverse metric access |
| `face_metrics[axis].forward` | face-centered forward Jacobian and determinant | diagnostics and face metric access |
| `face_metrics[axis].inverse` | face-centered inverse metric | diagnostics and face metric access |
| `face_metrics[axis].conserved` | conserved face metric matrix `Ghat` | solver-facing conservative face geometry |

Solvers should usually consume face geometry through:

```julia
geom = face_flux_geometry(grid, idx, :ihi)
geom.metric_vector
geom.area
geom.normal
```

`face_flux_geometry` reads the active row of
`face_metrics(grid)[axis].conserved` and applies coordinate-system measure
factors and outward orientation.

## Generic Reconstruction Path

For `EndpointAverageReconstruction`, `GradientCorrectedReconstruction`, and
`CurvatureCorrectedReconstruction`, the generic path is:

1. Build mapping derivative functions from the grid mapping.
2. Build inverse and conserved metric functions.
3. Reconstruct those metric functions to faces with the selected face
   reconstruction scheme.
4. Fill `face_metrics[axis].forward`, `face_metrics[axis].inverse`, and
   `face_metrics[axis].conserved`.

In 1D and 2D, the conserved face metric is essentially the reconstructed
`J * inv(F)` form. In 3D, the conserved metric functions follow the
Thomas-Lombard conservative form before reconstruction.

The selector controls the reconstruction used at the face:

| Selector | Face reconstruction data |
|---|---|
| `EndpointAverageReconstruction()` | endpoint values |
| `GradientCorrectedReconstruction()` | endpoint values and first derivatives |
| `CurvatureCorrectedReconstruction()` | endpoint values, first derivatives, and second derivatives |

For `DiscreteGrid`, this path remains important because the grid is defined by
sampled coordinates. Smooth global interpolation can oscillate, fold cells, or
invent geometry that was not present in the mesh. The conservative discrete
metric scheme works directly with the sampled mesh representation.

## AD-Driven Metrics

For `MappedGrid`, automatic differentiation can evaluate derivatives of the
actual analytic mapping. This is better than differentiating sampled coordinates
when the map is truly known.

`ADThomasLombardMetric()` combines AD-evaluated mapping information with
Thomas-Lombard conservative metric assembly in 3D. The intent is:

- preserve the conservative algebra needed for GCL,
- avoid finite-difference stencil error in local mapping derivatives,
- fill solver-facing conserved face metric rows consistently.

The important implementation detail is that `ADThomasLombardMetric()` is not
"compute conserved metrics at cell centers and interpolate them to faces."
Instead:

| Quantity | `ADThomasLombardMetric()` behavior |
|---|---|
| `cell_metrics.forward` and `cell_metrics.inverse` | AD Jacobian and inverse Jacobian evaluated at cell centers. This is the ordinary mapped-cell metric path. |
| `face_metrics[axis].forward` and `face_metrics[axis].inverse` | AD Jacobian and inverse Jacobian evaluated directly at face centers. |
| `face_metrics[axis].conserved` | Special Thomas-Lombard face assembly from AD-evaluated metric potentials. |
| stored conserved face matrix | Only the active row is populated: xi row on `i` faces, eta row on `j` faces, zeta row on `k` faces. |

The AD Thomas-Lombard conserved path evaluates metric potentials such as
products of coordinate derivatives and coordinates, reconstructs those potentials
at face locations, and accumulates paired derivatives so that the
Thomas-Lombard cancellation structure is retained.

In 1D and 2D, the AD Thomas-Lombard selector simplifies to the direct metric
form because the full 3D Thomas-Lombard pair structure is not needed. The
constructor emits a warning and uses the default direct scheme.

The special AD Thomas-Lombard face assembly is currently a CPU path. Treat it as
a 3D CPU `MappedGrid` option unless a backend-specific implementation or guard
is added.

## DiscreteGrid Metric Nuance

`DiscreteGrid` is interpolation-backed. Current constructors build linear
interpolants from node arrays. The metric cache wraps those interpolants as
mapping functions.

There is one important detail:

- forward cell Jacobians are overridden to use `Interpolations.gradient` on the
  linear interpolants,
- inverse and edge metric functions still follow the generic metric-cache path,
- AD or derivative calls on a `DiscreteGrid` differentiate the interpolation
  surrogate, not the original CAD or mesh-generation process.

For this reason, `ADThomasLombardMetric()` should not be treated as a magic
upgrade for imported 3D coordinates. On sampled meshes, the conservative
discrete reconstruction path is the safer default because the sampled nodes are
the geometry source of truth.

## Stencil Versus AD

| Question | Discrete stencil path | AD mapped path |
|---|---|---|
| Geometry source | sampled node coordinates | analytic mapping functions |
| Best use | CAD/imported/generated meshes | smooth exact maps |
| Differentiates | discrete coordinate representation | actual map |
| GCL strategy | conservative discrete operator compatibility | AD derivatives plus conservative assembly |
| Main risk | limited by mesh quality and interpolation smoothness | only valid when the map is genuinely known |
| Failure mode | bad/noisy/folded mesh still bad | smooth surrogate may not represent imported mesh |

The old rule of thumb was that the advection reconstruction and metric
derivation should use the same spatial reconstruction scheme. A more precise
modern statement is:

> The metric construction must be compatible with the discrete face locations,
> quadrature, and divergence operator used by the solver.

For analytic maps, AD can loosen the old stencil-matching rule because it
improves the derivative source. It does not remove the need for conservative
face assembly or GCL testing.

## Recommendations

Use these defaults unless a test or solver requirement says otherwise.

| Scenario | Recommended scheme | Rationale |
|---|---|---|
| 1D mapped or discrete grid | `CurvatureCorrectedReconstruction()` | Robust default. AD Thomas-Lombard falls back in 1D. |
| 2D mapped or discrete grid | `CurvatureCorrectedReconstruction()` | Robust default. AD Thomas-Lombard falls back in 2D. |
| 3D analytic `MappedGrid` on CPU | `ADThomasLombardMetric()` for metric-sensitive production tests; `CurvatureCorrectedReconstruction()` as reference/comparison | ADTL uses the analytic map and preserves Thomas-Lombard conservative assembly at faces. |
| 3D analytic `MappedGrid` on GPU | `CurvatureCorrectedReconstruction()` | The special AD Thomas-Lombard face assembly is currently CPU-specific. |
| CAD/imported/generated node coordinates | `DiscreteGrid(...; conserved_metric_scheme=CurvatureCorrectedReconstruction())` | Imported geometry is sampled data. Conservative discrete reconstruction is the trusted path. |
| Debugging a metric issue | Compare `EndpointAverageReconstruction()`, `GradientCorrectedReconstruction()`, and `CurvatureCorrectedReconstruction()` | Scheme differences help isolate reconstruction sensitivity. |
| Freestream preservation tests | Any candidate scheme must pass GCL and constant-state solver tests | GCL is necessary but not sufficient. |
| Riemann-solver accuracy on curved grids | Compare `face_flux_geometry` vectors, areas, normals, and coordinates | The Riemann solve consumes face geometry, not just cell metrics. |

Do not choose a scheme only because it has the highest local derivative accuracy.
The face metric must be compatible with the solver's face locations, divergence
operator, coordinate-system measure factors, and conservation update.

## What To Test When Changing Schemes

At minimum:

1. Check GCL residuals from `face_metrics(grid)`.
2. Check positive cell volumes and positive face areas.
3. Compare `face_flux_geometry` across refinement for known smooth geometries.
4. Run freestream preservation tests.
5. Run at least one non-constant smooth advection or manufactured-solution test
   on the target coordinate system.

For sampled/CAD grids, also check mesh-quality diagnostics: skewness, stretching,
normal smoothness, cell volume jumps, and surface curvature resolution.

## Current Mapped-Versus-Discrete Evidence

The test suite includes `MappedGrid` versus `DiscreteGrid` checks where the
mapped nodes are sampled into a discrete grid for the same geometry. The test
matrix covers valid 1D, 2D, and 3D coordinate-system pairings, checks GCL
residuals below `1e-12`, and verifies that solver-facing face geometry converges
under refinement.

The observed pattern is:

- 1D is essentially exact for the tested cases.
- 2D metric vectors and areas converge strongly, with normals and face
  coordinates converging more slowly.
- 3D metric vectors and areas still converge strongly, but normals and face
  coordinates are the weaker quantities.

That is the expected behavior for a `DiscreteGrid` backed by linear
interpolation: the conservative metric algebra is good, but the coordinate
surrogate is only piecewise smooth. This supports `DiscreteGrid` as the sampled
mesh path without claiming that it is AD-quality geometry at a fixed resolution.
