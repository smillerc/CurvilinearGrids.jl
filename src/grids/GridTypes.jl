module GridTypes

using LinearAlgebra
using StaticArrays
using ForwardDiff: derivative, gradient

export AbstractCurvilinearMesh
export CurvilinearMesh2D, CurvilinearMesh3D
export coord, coords
export centroid, centroids
export metrics, jacobian, jacobian_matrix

abstract type AbstractCurvilinearMesh end

include("2d.jl")
include("3d.jl")

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
