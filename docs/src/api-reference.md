# API Reference

This page is a compact index of the public API most relevant to unified grids
and solver integration.

## Grid Types

```julia
MappedGrid
DiscreteGrid
OrthogonalGrid

AbstractUnifiedGrid
Metric
ConservedMetric
FaceFluxGeometry
```

## Coordinate Systems and Bases

```julia
CoordinateSystemTrait
CurvilinearCS
CartesianCS
CylindricalCS
SphericalCS
AxisymmetricCS

BasisTrait
CartesianBasis
SphericalBasis

coordinate_system
basis_trait
basis_transfer_matrix
```

## Metric Schemes

```julia
EdgeInterpolationSchemeTrait
EndpointAverageReconstruction
GradientCorrectedReconstruction
CurvatureCorrectedReconstruction
ADThomasLombardMetric
```

## Geometry Access

```julia
coord
coords
cartesian_coordinates

centroid
centroids
cartesian_centroid
cartesian_centroids

face_coordinate
face_coordinates
face_flux_geometry
face_area
outward_face_normal

cellvolume
cellvolumes
cell_metrics
face_metrics
forward_cell_metrics
inverse_cell_metrics
jacobian_matrix
```

## Metric Cache Lifecycle

```julia
update!
invalidate_cell_metrics!
invalidate_face_metrics!
refresh_cell_metrics!
refresh_face_metrics!
```

## Coordinate-System Conservation Scales

```julia
conservation_cell_metric_scale
conservation_face_metric_component_scale
```

## IO and Utilities

```julia
save_vtk
write_coordinates
read_coordinates

computational_coordinate
InverseCoordinateResult
```
