module GridTypes

using LinearAlgebra
using StaticArrays
using ForwardDiff

export AbstractCurvilinearMesh
export CurvilinearMesh2D, CurvilinearMesh3D
export coord, coords
export centroid, centroids
export metrics, jacobian, jacobian_matrix

abstract type AbstractCurvilinearMesh end

# Helper functions to eliminate floating point values below ϵ; these
# are used to "clean" the jacobian matrices, so we don't get 
# numbers like 1e-18, when in reality they should be 0.0
# @inline checkeps(M) = M # non-Float version 
@inline function checkeps(M::AbstractArray{T,N}) where {T<:AbstractFloat,N}
  return @. M * abs(M >= eps(T))
end

# metric types, e.g. named tuples
Metrics2D{T} = NamedTuple{(:ξ̂x, :η̂x, :ζ̂x, :ξ̂y, :η̂y, :ζ̂y),NTuple{6,T}} where {T}
Metrics3D{T} =
  NamedTuple{(:ξ̂x, :η̂x, :ζ̂x, :ξ̂y, :η̂y, :ζ̂y, :ξ̂z, :η̂z, :ζ̂z),NTuple{9,T}} where {T}

function _make_empty_metric_array((ni, nj)::NTuple{2,Int}, T=Float64)
  M = Array{Metrics2D{Float64},2}(undef, ni, nj)

  @inbounds for i in eachindex(M)
    M[i] = (ξ̂x=zero(T), η̂x=zero(T), ζ̂x=zero(T), ξ̂y=zero(T), η̂y=zero(T), ζ̂y=zero(T))
  end

  return M
end

function _make_empty_metric_array((ni, nj, nk)::NTuple{3,Int}, T=Float64)
  M = Array{Metrics3D{Float64},3}(undef, ni, nj, nk)

  @inbounds for i in eachindex(M)
    M[i] = (
      ξ̂x=zero(T),
      η̂x=zero(T),
      ζ̂x=zero(T),
      ξ̂y=zero(T),
      η̂y=zero(T),
      ζ̂y=zero(T),
      ξ̂z=zero(T),
      η̂z=zero(T),
      ζ̂z=zero(T),
    )
  end

  return M
end

include("2d.jl")
include("3d.jl")

function update(m::AbstractCurvilinearMesh)
  #TODO: make these GPU/CPU agnostic
  update_metrics(m)

  @inbounds for I in CartesianIndices(m.J)
    _jacobian_matrix = jacobian_matrix(m, I)
    m.J[I] = det(_jacobian_matrix)
    # m.J⁻¹[I] = det(inv(_jacobian_matrix))
  end

  return nothing
end

coord(mesh, CI::CartesianIndex) = coord(mesh, CI.I...)
coord(m::CurvilinearMesh2D, i, j) = @SVector [m.x(i, j), m.y(i, j)]
coord(m::CurvilinearMesh3D, i, j, k) = @SVector [m.x(i, j, k), m.y(i, j, k), m.z(i, j, k)]

centroid(mesh, CI::CartesianIndex) = centroid(mesh, CI.I...)

function centroid(m::CurvilinearMesh2D, i, j)
  @SVector [m.x(i + 0.5, j + 0.5), m.y(i + 0.5, j + 0.5)]
end

function centroid(m::CurvilinearMesh3D, i, j, k)
  @SVector [
    m.x(i + 0.5, j + 0.5, k + 0.5),
    m.y(i + 0.5, j + 0.5, k + 0.5),
    m.z(i + 0.5, j + 0.5, k + 0.5),
  ]
end

jacobian_matrix(mesh, CI::CartesianIndex) = mesh.jacobian_matrix_func(CI.I...)
jacobian_matrix(mesh, idx) = jacobian_matrix(mesh, idx...)
jacobian_matrix(mesh, i, j, k) = mesh.jacobian_matrix_func(i, j, k)
jacobian_matrix(mesh, i, j) = mesh.jacobian_matrix_func(i, j)

jacobian(mesh, CI::CartesianIndex) = det(jacobian_matrix(mesh, CI))
jacobian(mesh, idx) = det(jacobian_matrix(mesh, idx))
jacobian(mesh::CurvilinearMesh3D, i, j, k) = det(jacobian_matrix(mesh, i, j, k))
jacobian(mesh::CurvilinearMesh2D, i, j) = det(jacobian_matrix(mesh, i, j))

"""
Query the mesh metrics at a particular index

**Note**: The mesh stores these in `mesh.cell_center_metrics` and `mesh.edge_metrics`

"""
metrics(mesh, CI::CartesianIndex) = metrics(mesh, CI.I...)
metrics(mesh, CI::CartesianIndex, v⃗_grid) = metrics(mesh, CI.I..., v⃗_grid)

metrics(mesh, idx) = metrics(mesh, idx...)
metrics(mesh, idx, v⃗_grid) = metrics(mesh, idx..., v⃗_grid)

function metrics(m::CurvilinearMesh3D, i, j, k)
  jacobi = m.jacobian_matrix_func(i, j, k)
  return _get_metrics(jacobi)
end

function metrics(m::CurvilinearMesh2D, i, j)
  jacobi = m.jacobian_matrix_func(i, j)
  return _get_metrics(jacobi)
end

function metrics(m::CurvilinearMesh3D, i, j, k, v⃗_grid)
  jacobi = m.jacobian_matrix_func(i, j, k)
  return _get_metrics(jacobi, v⃗_grid)
end

function metrics(m::CurvilinearMesh2D, i, j, v⃗_grid)
  jacobi = m.jacobian_matrix_func(i, j)
  return _get_metrics(jacobi, v⃗_grid)
end

end
