module CurvilinearGrids

include("grids/GridTypes.jl")
using .GridTypes
export CurvilinearMesh2D
export coords, xy, centroid_xy, centroids, metrics

include("io")

end
