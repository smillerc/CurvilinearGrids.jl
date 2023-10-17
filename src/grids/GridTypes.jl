module GridTypes

using LinearAlgebra
using StaticArrays
using ForwardDiff
using UnPack

export AbstractCurvilinearGrid
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export coord, coords, cellsize
export centroid, centroids
export metrics, jacobian, jacobian_matrix
export conservative_metrics

abstract type AbstractCurvilinearGrid end

# Helper functions to eliminate floating point values below ϵ; these
# are used to "clean" the jacobian matrices, so we don't get 
# numbers like 1e-18, when in reality they should be 0.0
# @inline checkeps(M) = M # non-Float version 
@inline function checkeps(M::AbstractArray{T}) where {T<:AbstractFloat}
  return @. M * abs(M >= eps(T))
end
@inline checkeps(M) = M

include("1d.jl")
include("2d.jl")
include("3d.jl")

"""Get the size of the grid for cell-based arrays"""
cellsize(m::AbstractCurvilinearGrid) = @. m.nnodes - 1

"""Get the size of the grid for cell-based arrays when the halo cells are included"""
cellsize_withhalo(m::AbstractCurvilinearGrid) = @. m.nnodes - 1 + 2 * m.nhalo

"""

Get the position/coordinate at a given index. **NOTE:** these indices are consistent with halo cells included.
This means that if your grid has 2 halo cells, the position of the first non-halo vertex is 
index at coord(mesh, 3). The `CurvilinearGrid` only keeps track of the number of halo cells for each
dimension, whereas the grid functions have no knowledge halos. Therefore, the `coord` function
applies a shift to the index for you.
"""
coord(mesh, CI::CartesianIndex) = coord(mesh, CI.I...)
coord(m::CurvilinearGrid1D, i) = m.x(i)
coord(m::CurvilinearGrid2D, i, j) = @SVector [m.x(i, j), m.y(i, j)]
coord(m::CurvilinearGrid3D, i, j, k) = @SVector [m.x(i, j, k), m.y(i, j, k), m.z(i, j, k)]

"""

Get the position of the centroid for the given _cell_ index. **NOTE:** these indices 
are consistent with halo cells included. This means that if your grid has 2 halo cells, 
the position of the first non-halo centroid is index at coord(mesh, 3). 
The `CurvilinearGrid` only keeps track of the number of halo cells for each dimension, 
whereas the grid functions have no knowledge halos. Therefore, the `coord` function
applies a shift to the index for you.
"""
centroid(mesh, CI::CartesianIndex) = centroid(mesh, CI.I...)
centroid(m::CurvilinearGrid1D, i) = m.x(i - m.nhalo + 0.5)

function centroid(m::CurvilinearGrid2D, i, j)
  @SVector [
    m.x(i - m.nhalo + 0.5, j - m.nhalo + 0.5), # x
    m.y(i - m.nhalo + 0.5, j - m.nhalo + 0.5), # y
  ]
end

function centroid(m::CurvilinearGrid3D, i, j, k)
  @SVector [
    m.x(i - m.nhalo + 0.5, j - m.nhalo + 0.5, k - m.nhalo + 0.5), # x
    m.y(i - m.nhalo + 0.5, j - m.nhalo + 0.5, k - m.nhalo + 0.5), # y
    m.z(i - m.nhalo + 0.5, j - m.nhalo + 0.5, k - m.nhalo + 0.5), # z
  ]
end

"""
Get the Jacobian matrix of the forward transformation (ξ,η,ζ) → (x,y,z).
"""
jacobian_matrix(mesh, CI::CartesianIndex) = jacobian_matrix(mesh, CI.I)
jacobian_matrix(mesh, idx::NTuple) = jacobian_matrix(mesh, idx...)
function jacobian_matrix(mesh::CurvilinearGrid3D, i, j, k)
  return checkeps(mesh.jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo, k - mesh.nhalo))
end

function jacobian_matrix(mesh::CurvilinearGrid2D, i, j)
  return checkeps(mesh.jacobian_matrix_func(i - mesh.nhalo, j - mesh.nhalo))
end

"""
Get the Jacobian of the forward transformation (ξ,η,ζ) → (x,y,z).
"""
jacobian(mesh, CI::CartesianIndex) = det(jacobian_matrix(mesh, CI.I))
jacobian(mesh, idx::NTuple) = jacobian(mesh, idx...)
jacobian(mesh::CurvilinearGrid3D, i, j, k) = det(jacobian_matrix(mesh, (i, j, k)))
jacobian(mesh::CurvilinearGrid2D, i, j) = det(jacobian_matrix(mesh, (i, j)))
jacobian(mesh::CurvilinearGrid1D, i) = det(jacobian_matrix(mesh, i))

"""
Query the mesh metrics at a particular index
"""
metrics(mesh, CI::CartesianIndex) = metrics(mesh, CI.I...)
metrics(mesh, idx::NTuple{N,T}) where {N,T} = metrics(mesh, idx...)

end
