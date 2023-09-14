module GridTypes

using LinearAlgebra
using StaticArrays
using ForwardDiff: derivative, gradient

export CurvilinearMesh2D
export coords, xy, centroid_xy, centroids, metrics

abstract type AbstractCurvilinearMesh end

include("2d.jl")
end