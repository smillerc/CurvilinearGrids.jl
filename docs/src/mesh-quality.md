# Mesh Quality and Testing

Metric schemes can preserve conservation, but they cannot turn a bad mesh into a
good one. This is especially important for `DiscreteGrid`, which is the current
path for CAD-generated and imported coordinates.

## Practical Mesh Checks

Any production CFD mesh should be checked for:

| Category | Check |
|---|---|
| validity | positive cell volumes/Jacobians and positive face areas |
| smoothness | bounded cell-to-cell changes in spacing, area, volume, and normal |
| skewness | controlled non-orthogonality and face-center offset |
| aspect ratio | expected stretching without abrupt jumps |
| curvature | adequate surface and normal resolution |
| conservation | GCL residuals near roundoff for mapped/discrete metrics |
| solver geometry | sane `face_flux_geometry` vectors, areas, normals, and coordinates |

## Mapped Versus Discrete Regression Tests

A useful baseline test is:

1. Build a smooth `MappedGrid`.
2. Sample its node coordinates.
3. Build a `DiscreteGrid` from those nodes.
4. Check GCL residuals for both grids.
5. Compare `face_flux_geometry` outputs.
6. Repeat under refinement and require convergence.

This proves more than freestream preservation. It checks the actual face vectors,
areas, normals, and face coordinates that solvers consume.

The CurvilinearGrids test suite includes this matrix for valid 1D, 2D, and 3D
coordinate-system cases. The current behavior is:

| Quantity | Assessment |
|---|---|
| GCL | near roundoff for both mapped and discrete grids |
| face metric vector | strong convergence on smooth sampled grids |
| face area | strong convergence on smooth sampled grids |
| face normal | converges more slowly, especially in 3D |
| face coordinate | weakest convergence for the linear discrete surrogate |

This is a good result for sampled coordinates. It means the conservative
discrete path is viable, but it should not be mistaken for an exact AD
representation of the original smooth map.

## CAD and Imported Meshes

For CAD/imported coordinates, `DiscreteGrid` is the appropriate route because
the mesh usually provides nodes, not an analytic differentiable map. Fitting a
smooth multidimensional interpolant just to use AD can introduce oscillations,
folded cells, or nonlocal artifacts. The safer default is:

- keep the sampled coordinates as the source of truth,
- use conservative metric reconstruction,
- validate mesh quality explicitly,
- compare face areas/normals against independent geometric references when
  possible.

AD remains the best path when the analytic mapping is actually known.
