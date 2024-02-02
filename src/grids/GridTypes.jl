module GridTypes

using LinearAlgebra
using StaticArrays
using ForwardDiff
using UnPack

using ..MetricTypes
using ..MetricDiscretizationSchemes

export AbstractCurvilinearGrid
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export CylindricalGrid1D, SphericalGrid1D
export RZAxisymmetricGrid2D
export coord, coords, coords!, cellsize, cellsize_withhalo
export centroid, centroids
export area, volume
export metrics, jacobian, jacobian_matrix
export conservative_metrics
export metrics_with_jacobian
export cell_metrics, cell_indices

abstract type AbstractCurvilinearGrid end

# Helper functions to eliminate floating point values below ϵ; these
# are used to "clean" the jacobian matrices, so we don't get
# numbers like 1e-18, when in reality they should be 0.0
# @inline checkeps(M) = M # non-Float version
@inline function checkeps(M::AbstractArray{T}) where {T<:AbstractFloat}
  # return M #@. M * (abs(M) >= eps(T))
  return @. M * (abs(M) >= eps(T))
end
@inline checkeps(M) = M

# helper func to ensure that coordinate functions
# take the proper number of arguments, e.g., x(i) for 1d, x(i,j) for 2d
function check_nargs(f, nargs, fname)
  for m in methods(f)
    if (m.nargs - 1) != nargs
      error(
        "The function $(fname)() isn't defined correctly, it needs to be use $nargs (not $(m.nargs - 1)) arguments, e.g. `x(i,j) = ...` for 2D",
      )
    end
  end
end

# helper func to ensure that coordinate functions are working properly
function test_coord_func(f, ndims, fname)
  args = ones(Int, ndims)
  try
    f(args...)
  catch err
    @error(
      "Unable to evaluate the coordinate function $fname, due to the following errors:"
    )
    throw(err)
  end
end

include("1d.jl")
include("1d_axisymmetric.jl")
include("2d.jl")
include("2d_axisymmetric.jl")
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
coord(mesh, CI::CartesianIndex) = coord(mesh, CI.I)
coord(m::CurvilinearGrid1D, i) = m.x(i - m.nhalo)
coord(m::CurvilinearGrid1D, (i,)::NTuple{1,Real}) = m.x(i - m.nhalo)

function coord(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  @SVector [m.x(i - m.nhalo, j - m.nhalo), m.y(i - m.nhalo, j - m.nhalo)]
end

function coord(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})
  @SVector [m.r(i - m.nhalo, 1, j - m.nhalo), m.z(i - m.nhalo, 1, j - m.nhalo)]
end

function coord(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  @SVector [
    m.r(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m.θ(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m.z(i - m.nhalo, j - m.nhalo, k - m.nhalo),
  ]
end

function coord(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  @SVector [
    m.x(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m.y(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m.z(i - m.nhalo, j - m.nhalo, k - m.nhalo),
  ]
end

"""

Get the position of the centroid for the given _cell_ index. **NOTE:** these indices
are consistent with halo cells included. This means that if your grid has 2 halo cells,
the position of the first non-halo centroid is index at coord(mesh, 3).
The `CurvilinearGrid` only keeps track of the number of halo cells for each dimension,
whereas the grid functions have no knowledge halos. Therefore, the `coord` function
applies a shift to the index for you.
"""
centroid(mesh, CI::CartesianIndex) = centroid(mesh, CI.I)
centroid(m::CurvilinearGrid1D, i) = m.x(i - m.nhalo + 0.5)

function centroid(m::CurvilinearGrid2D, (i, j)::NTuple{2,Int})
  @SVector [
    m.x_func(i - m.nhalo + 0.5, j - m.nhalo + 0.5), # x
    m.y_func(i - m.nhalo + 0.5, j - m.nhalo + 0.5), # y
  ]
end

function centroid(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Int})
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
jacobian_matrix(mesh::CurvilinearGrid1D, i) = jacobian_matrix(mesh, (i,))

"""
Get the Jacobian of the forward transformation (ξ,η,ζ) → (x,y,z).
"""
jacobian(mesh, CI::CartesianIndex) = jacobian(mesh, CI.I)
jacobian(mesh::CurvilinearGrid1D, i::Real) = jacobian(mesh, (i,))

"""
Query the mesh metrics at a particular index
"""
metrics(mesh, CI::CartesianIndex) = metrics(mesh, CI.I...)
metrics(mesh::CurvilinearGrid1D, i::Real) = metrics(mesh, (i,))

@inline function cell_metrics(m::CurvilinearGrid1D, i::Real)
  return metrics(m, i + 0.5)
end

@inline function cell_metrics(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  return metrics(m, (i + 0.5, j + 0.5))
end

@inline function cell_metrics(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  return metrics(m, (i + 0.5, j + 0.5, k + 0.5))
end

"""
Query the conservative mesh metrics at a particular index that follow the GCL
"""
conservative_metrics(mesh, CI::CartesianIndex) = conservative_metrics(mesh, CI.I...)
conservative_metrics(mesh::CurvilinearGrid1D, i::Real) = conservative_metrics(mesh, (i,))

function cell_indices(mesh::CurvilinearGrid1D)
  @unpack ilo, ihi = mesh.limits
  return CartesianIndices((ilo:ihi))
end

function cell_indices(mesh::CurvilinearGrid2D)
  @unpack ilo, ihi, jlo, jhi = mesh.limits
  return CartesianIndices((ilo:ihi, jlo:jhi))
end

function cell_indices(mesh::CurvilinearGrid3D)
  @unpack ilo, ihi, jlo, jhi, klo, khi = mesh.limits
  return CartesianIndices((ilo:ihi, jlo:jhi, klo:khi))
end

end
