module GridTypes

using LinearAlgebra
using StaticArrays
using ForwardDiff
using UnPack
using StructArrays
using KernelAbstractions

using ..IndexingUtils
using ..MetricTypes
using ..MetricDiscretizationSchemes

export AbstractCurvilinearGrid
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export CylindricalGrid1D, SphericalGrid1D
export RZAxisymmetricGrid2D
export coord, coords, coords!, cellsize, cellsize_withhalo
export centroid, centroids
export cellvolume
export metrics, jacobian, jacobian_matrix
export conservative_metrics
export metrics_with_jacobian
export cell_metrics, cell_indices

abstract type AbstractCurvilinearGrid end

# Helper functions to eliminate floating point values below ϵ; these
# are used to "clean" the jacobian matrices, so we don't get
# numbers like 1e-18, when in reality they should be 0.0
# @inline checkeps(M) = M # non-Float version
@inline function checkeps(M::AbstractArray{T}, ϵ=eps(T)) where {T<:AbstractFloat}
  # return M #@. M * (abs(M) >= eps(T))
  return @. M * (abs(M) >= ϵ)
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

include("metric_soa.jl")
include("grid_iterators.jl")
include("1d.jl")
# include("1d_spherical.jl")
include("1d_axisymmetric.jl")
include("2d.jl")
include("2d_axisymmetric.jl")
include("3d.jl")

"""Get the size of the grid for cell-based arrays"""
cellsize(m::AbstractCurvilinearGrid) = @. m.nnodes - 1

"""Get the size of the grid for cell-based arrays when the halo cells are included"""
cellsize_withhalo(m::AbstractCurvilinearGrid) = @. m.nnodes - 1 + 2 * m.nhalo

function coords(m::CurvilinearGrid1D)
  return @views m.node_coordinates.x[m.iterators.node.domain]
end

function coords(m::CurvilinearGrid2D)
  return @views (
    m.node_coordinates.x[m.iterators.node.domain],
    m.node_coordinates.y[m.iterators.node.domain],
  )
end

function coords(m::CurvilinearGrid3D)
  return @views (
    m.node_coordinates.x[m.iterators.node.domain],
    m.node_coordinates.y[m.iterators.node.domain],
    m.node_coordinates.z[m.iterators.node.domain],
  )
end

function centroids(m::CurvilinearGrid1D)
  return @views m.centroid_coordinates.x[m.iterators.cell.domain]
end

function centroids(m::CurvilinearGrid2D)
  return @views (
    m.centroid_coordinates.x[m.iterators.cell.domain],
    m.centroid_coordinates.y[m.iterators.cell.domain],
  )
end

function centroids(m::CurvilinearGrid3D)
  return @views (
    m.centroid_coordinates.x[m.iterators.cell.domain],
    m.centroid_coordinates.y[m.iterators.cell.domain],
    m.centroid_coordinates.z[m.iterators.cell.domain],
  )
end

"""

Get the position/coordinate at a given index. **NOTE:** these indices are consistent with halo cells included.
This means that if your grid has 2 halo cells, the position of the first non-halo vertex is
index at coord(mesh, 3). The `CurvilinearGrid` only keeps track of the number of halo cells for each
dimension, whereas the grid functions have no knowledge halos. Therefore, the `coord` function
applies a shift to the index for you.
"""
coord(mesh, CI::CartesianIndex) = coord(mesh, CI.I)
coord(m, i::Real) = coord(m, (i,))
coord(m::CurvilinearGrid1D, (i,)::NTuple{1,Real}) = m._coordinate_funcs.x(i - m.nhalo)
coord(m::AxisymmetricGrid1D, (i,)::NTuple{1,Real}) = m._coordinate_funcs.r(i - m.nhalo)

function coord(m::CurvilinearGrid2D, (i, j)::NTuple{2,Real})
  @SVector [
    m._coordinate_funcs.x(i - m.nhalo, j - m.nhalo),
    m._coordinate_funcs.y(i - m.nhalo, j - m.nhalo),
  ]
end

function coord(m::RZAxisymmetricGrid2D, (i, j)::NTuple{2,Real})
  @SVector [
    m._coordinate_funcs.r(i - m.nhalo, 1, j - m.nhalo),
    m._coordinate_funcs.z(i - m.nhalo, 1, j - m.nhalo),
  ]
end

function coord(m::RZAxisymmetricGrid2D, (i, j, k)::NTuple{3,Real})
  @SVector [
    m._coordinate_funcs.r(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m._coordinate_funcs.θ(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m._coordinate_funcs.z(i - m.nhalo, j - m.nhalo, k - m.nhalo),
  ]
end

function coord(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real})
  @SVector [
    m._coordinate_funcs.x(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m._coordinate_funcs.y(i - m.nhalo, j - m.nhalo, k - m.nhalo),
    m._coordinate_funcs.z(i - m.nhalo, j - m.nhalo, k - m.nhalo),
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
centroid(m, i::Real) = coord(m, i + 0.5)

function centroid(m::CurvilinearGrid2D, (i, j)::NTuple{2,Int})
  return coord(m, (i, j) .+ 0.5)
end

function centroid(m::CurvilinearGrid3D, (i, j, k)::NTuple{3,Int})
  return coord(m, (i, j, k) .+ 0.5)
end

"""
Get the Jacobian matrix of the forward transformation (ξ,η,ζ) → (x,y,z).
"""
jacobian_matrix(mesh, CI::CartesianIndex, t::Real=0) = jacobian_matrix(mesh, CI.I, t)
jacobian_matrix(mesh, i, t::Real=0) = jacobian_matrix(mesh, (i,), t)

"""
Get the Jacobian of the forward transformation (ξ,η,ζ) → (x,y,z).
"""
jacobian(mesh, CI::CartesianIndex, t::Real=0) = jacobian(mesh, CI.I, t)
jacobian(mesh, i::Real, t::Real=0) = jacobian(mesh, (i,), t)

"""
Query the mesh metrics at a particular index
"""
@inline metrics(mesh, CI::CartesianIndex, t::Real=0) = metrics(mesh, CI.I, t)
@inline metrics(mesh, i::Real, t::Real=0) = metrics(mesh, (i,), t)

"""
Query the conservative mesh metrics at a particular index that follow the GCL
"""
@inline conservative_metrics(mesh, CI::CartesianIndex) = conservative_metrics(mesh, CI.I, t)
@inline conservative_metrics(mesh, i::Real, t::Real=0) = conservative_metrics(mesh, (i,), t)

function check_for_invalid_metrics(m::AbstractCurvilinearGrid)
  domain = m.iterators.cell.domain

  # cell centroid metrics
  invalid_cell_metrics = false
  for idx in domain
    if !(isfinite(m.cell_center_metrics.J[idx]))
      @error("Invalid jacobian @ index $(idx)")
      invalid_cell_metrics = true
    end
  end

  for (name, data) in pairs(m.cell_center_metrics)
    if data isa StructArray
      for idx in domain
        for c in StructArrays.components(data)
          if !(isfinite(c[idx]))
            @error("Invalid grid metric $(name) @ index $(idx)")
            invalid_cell_metrics = true
          end
        end
      end
    end
  end

  # The edge metrics have a domain expanded by 1 along
  # the particular dimension, since they store the edge "i+1/2"
  # for each cell. The lowest-most cell needs an i+1/2 too, so the
  # cell at ilo - 1 needs to have a valid value, where ilo is the first
  # index of the cell in the non-halo domain along the i axis.
  # iaxis, jaxis, kaxis = (1, 2, 3)
  # i₊½CI = expand_lower(domain, iaxis, 1)
  # j₊½CI = expand_lower(domain, jaxis, 1)
  # k₊½CI = expand_lower(domain, kaxis, 1)
  edge_iterators = ntuple(i -> expand_lower(domain, i, 1), length(m.nnodes))

  # edge metrics
  invalid_edge_metrics = false
  for (edge, edge_indices) in zip(m.edge_metrics, edge_iterators)
    for (name, data) in pairs(edge) # i₊½, j₊½, k₊½
      if data isa StructArray
        # metric_set = edge[idx] # ξ̂x, η̂x, etc...
        for c in StructArrays.components(data)
          for i in eachindex(edge_indices)
            if !(isfinite(c[i]))
              @error("Invalid conserved grid metric $(name) @ index $(i) of $(c[i])")
              invalid_edge_metrics = true
            end
          end
        end
      else
        for i in eachindex(edge_indices)
          if !(isfinite(data[i]))
            @error("Invalid conserved grid metric $(name) @ index $(i) of $(data[i])")
            invalid_edge_metrics = true
          end
        end
      end
    end
  end

  if invalid_cell_metrics || invalid_edge_metrics
    error(
      "Invalide cell or edge grid metrics, see error messages for what the culprits are; Exiting...",
    )
  end

  return nothing
end

end
