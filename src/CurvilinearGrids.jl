module CurvilinearGrids

include("grids/GridTypes.jl")
using .GridTypes
export CurvilinearMesh2D, CurvilinearMesh3D
export coord, coords
export centroid, centroids
export metrics, jacobian, jacobian_matrix

include("io/to_vtk.jl")
using .VTKOutput
export to_vtk

end
