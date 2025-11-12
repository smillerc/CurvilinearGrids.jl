module CurvilinearGrids

include("discretization_schemes/DiscretizationSchemes.jl")
using .DiscretizationSchemes
export MonotoneExplicitGradientScheme
export SecondOrder, FourthOrder, SixthOrder, EighthOrder
export compute_first_derivatives!, compute_second_derivatives!
export cell_center_derivatives!, interpolate_to_edge!

include("metric_schemes/MetricDiscretizationSchemes.jl")
using .MetricDiscretizationSchemes

include("RectilinearArrays.jl")
using .RectilinearArrays

include("grids/GridTypes.jl")
using .GridTypes
export AbstractCurvilinearGrid
export AbstractCurvilinearGrid1D
export AbstractCurvilinearGrid2D
export AbstractCurvilinearGrid3D
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export RectilinearGrid2D, RectilinearGrid3D
export UniformGrid1D, UniformGrid2D, UniformGrid3D
export CylindricalGrid1D, SphericalGrid1D
export AxisymmetricGrid2D

export ContinuousCurvilinearGrid1D, ContinuousCurvilinearGrid2D, ContinuousCurvilinearGrid3D

export rectilinear_grid,
  rtheta_grid, rthetaphi_grid, rectilinear_cylindrical_grid, rectilinear_spherical_grid
export axisymmetric_rectilinear_grid, axisymmetric_rtheta_grid

export update!
export cellsize, cellsize_withhalo
export coord, coords
export centroid, centroids
export metrics, jacobian, jacobian_matrix
export conservative_metrics
export metrics_with_jacobian
export cell_metrics, cell_indices
export radius, centroid_radius, centroid_radii
export cellvolume, cellvolumes

export forward_cell_metrics, inverse_cell_metrics

include("adapt_to_gpu.jl")

include("io/to_vtk.jl")
using .VTKOutput
export save_vtk

include("mesh_functions/stretching_functions.jl")
export one_sided_stretch, double_sided_stretch, one_sided_with_initial_spacing

include("grids/surface_mesh.jl")
export extract_surface_mesh

include("mesh_functions/estimate_yplus.jl")
export estimate_wall_distance

include("mesh_functions/rigid_body_transformations.jl")
export translate!, rotate!, scale!

include("remap.jl")
export change_resolution, scale_resolution, remap_cell_data

end
