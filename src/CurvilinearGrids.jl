module CurvilinearGrids

include("grids/GridTypes.jl")
using .GridTypes
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export cellsize
export coord, coords
export centroid, centroids
export metrics, jacobian, jacobian_matrix
export conservative_metrics

end
