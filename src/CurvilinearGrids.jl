module CurvilinearGrids

include("metrics.jl")
using .MetricTypes

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
export RectlinearGrid, RThetaGrid, RThetaPhiGrid
export AxisymmetricRectlinearGrid, AxisymmetricRThetaGrid

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

include("adapt_to_gpu.jl")

include("io/to_vtk.jl")
using .VTKOutput
export save_vtk

end
