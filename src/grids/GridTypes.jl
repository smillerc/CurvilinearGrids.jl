module GridTypes

using LinearAlgebra
using StaticArrays
using ForwardDiff: derivative, gradient

export AbstractCurvilinearMesh
export CurvilinearMesh2D, CurvilinearMesh3D
export coords, xy, centroid_xy, centroids, metrics

abstract type AbstractCurvilinearMesh end

include("2d.jl")
include("3d.jl")

end
