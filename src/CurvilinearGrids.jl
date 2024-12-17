module CurvilinearGrids

include("discretization_schemes/MetricDiscretizationSchemes.jl")
using .MetricDiscretizationSchemes

include("grids/GridTypes.jl")
using .GridTypes
export AbstractCurvilinearGrid
export AbstractCurvilinearGrid1D
export AbstractCurvilinearGrid2D
export AbstractCurvilinearGrid3D
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export CylindricalGrid1D, SphericalGrid1D
export AxisymmetricGrid2D
export rectlinear_grid,
  rtheta_grid, rthetaphi_grid, rectlinear_cylindrical_grid, rectlinear_spherical_grid
export axisymmetric_rectlinear_grid, axisymmetric_rtheta_grid

export update!
export cellsize, cellsize_withhalo
export coord, coords
export centroid, centroids
export metrics, jacobian, jacobian_matrix
export conservative_metrics
export metrics_with_jacobian
export cell_metrics, cell_indices
export radius, centroid_radius
export cellvolume, cellvolumes

export forward_cell_metrics, inverse_cell_metrics

include("adapt_to_gpu.jl")

include("io/to_vtk.jl")
using .VTKOutput
export save_vtk

end
