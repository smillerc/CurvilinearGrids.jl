module CurvilinearGrids

include("grids/GridTypes.jl")
using .GridTypes
export CurvilinearMesh2D, CurvilinearMesh3D
export coords, xy, centroid_xy, centroids, metrics

include("io/to_vtk.jl")
using .VTKOutput
export to_vtk

end
