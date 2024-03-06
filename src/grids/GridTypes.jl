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
export CylindricalGrid2D
export coord, coords, coords!, cellsize, cellsize_withhalo
export centroid, centroids
export cellvolume
export radius, centroid_radius
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
@inline function check_nargs(f, nargs, fname)
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
cellsize(mesh::AbstractCurvilinearGrid) = @. mesh.nnodes - 1

"""Get the size of the grid for cell-based arrays when the halo cells are included"""
cellsize_withhalo(mesh::AbstractCurvilinearGrid) = @. mesh.nnodes - 1 + 2 * mesh.nhalo

@inline function coords(mesh::CurvilinearGrid1D)
  return @views mesh.node_coordinates.x[mesh.iterators.node.domain]
end

@inline function coords(mesh::CurvilinearGrid2D)
  return @views (
    mesh.node_coordinates.x[mesh.iterators.node.domain],
    mesh.node_coordinates.y[mesh.iterators.node.domain],
  )
end

@inline function coords(mesh::CylindricalGrid2D)
  return @views (
    mesh.node_coordinates.r[mesh.iterators.node.domain],
    mesh.node_coordinates.z[mesh.iterators.node.domain],
  )
end

@inline function coords(mesh::CurvilinearGrid3D)
  return @views (
    mesh.node_coordinates.x[mesh.iterators.node.domain],
    mesh.node_coordinates.y[mesh.iterators.node.domain],
    mesh.node_coordinates.z[mesh.iterators.node.domain],
  )
end

@inline function centroids(mesh::CurvilinearGrid1D)
  return @views mesh.centroid_coordinates.x[mesh.iterators.cell.domain]
end

@inline function centroids(mesh::CurvilinearGrid2D)
  return @views (
    mesh.centroid_coordinates.x[mesh.iterators.cell.domain],
    mesh.centroid_coordinates.y[mesh.iterators.cell.domain],
  )
end

@inline function centroids(mesh::CylindricalGrid2D)
  return @views (
    mesh.centroid_coordinates.r[mesh.iterators.cell.domain],
    mesh.centroid_coordinates.z[mesh.iterators.cell.domain],
  )
end

@inline function centroids(mesh::CurvilinearGrid3D)
  return @views (
    mesh.centroid_coordinates.x[mesh.iterators.cell.domain],
    mesh.centroid_coordinates.y[mesh.iterators.cell.domain],
    mesh.centroid_coordinates.z[mesh.iterators.cell.domain],
  )
end

"""

Get the position/coordinate at a given index. **NOTE:** these indices are consistent with halo cells included.
This means that if your grid has 2 halo cells, the position of the first non-halo vertex is
index at coord(mesh, 3). The `CurvilinearGrid` only keeps track of the number of halo cells for each
dimension, whereas the grid functions have no knowledge halos. Therefore, the `coord` function
applies a shift to the index for you.
"""
@inline coord(mesh, CI::CartesianIndex, t::Real=0) = coord(mesh, CI.I, t)

@inline function coord(mesh::CurvilinearGrid1D, (i,)::NTuple{1,Real}, t::Real=0)
  return mesh._coordinate_funcs.x(i - mesh.nhalo)
end

@inline function coord(mesh::AxisymmetricGrid1D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  return mesh._coordinate_funcs.r(i - mesh.nhalo)
end

@inline function coord(mesh::AxisymmetricGrid1D, (i,)::NTuple{1,Real}, t::Real=0)
  return mesh._coordinate_funcs.r(i - mesh.nhalo)
end

@inline function coord(mesh::CurvilinearGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  @SVector [
    mesh._coordinate_funcs.x(i - mesh.nhalo, j - mesh.nhalo),
    mesh._coordinate_funcs.y(i - mesh.nhalo, j - mesh.nhalo),
  ]
end

@inline function coord(mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  @SVector [
    mesh._coordinate_funcs.r(i - mesh.nhalo, j - mesh.nhalo),
    mesh._coordinate_funcs.z(i - mesh.nhalo, j - mesh.nhalo),
  ]
end

# @inline function coord(mesh::CylindricalGrid2D, (i, j, k)::NTuple{3,Real}, t::Real=0)
#   @SVector [
#     mesh._coordinate_funcs.r(i - mesh.nhalo, j - mesh.nhalo, k - mesh.nhalo),
#     mesh._coordinate_funcs.θ(i - mesh.nhalo, j - mesh.nhalo, k - mesh.nhalo),
#     mesh._coordinate_funcs.z(i - mesh.nhalo, j - mesh.nhalo, k - mesh.nhalo),
#   ]
# end

@inline function coord(mesh::CurvilinearGrid3D, (i, j, k)::NTuple{3,Real}, t::Real=0)
  @SVector [
    mesh._coordinate_funcs.x(i - mesh.nhalo, j - mesh.nhalo, k - mesh.nhalo),
    mesh._coordinate_funcs.y(i - mesh.nhalo, j - mesh.nhalo, k - mesh.nhalo),
    mesh._coordinate_funcs.z(i - mesh.nhalo, j - mesh.nhalo, k - mesh.nhalo),
  ]
end

@inline centroid_radius(mesh, CI::CartesianIndex, t::Real=0) = radius(mesh, CI.I .+ 0.5, t)
@inline radius(mesh, CI::CartesianIndex, t::Real=0) = radius(mesh, CI.I, t)

@inline function radius(mesh::AxisymmetricGrid1D, (i,)::NTuple{1,Real}, t::Real=0)
  return mesh._coordinate_funcs.r(i - mesh.nhalo) |> abs
end

@inline function radius(mesh::CylindricalGrid2D, (i, j)::NTuple{2,Real}, t::Real=0)
  return mesh._coordinate_funcs.r(i - mesh.nhalo, j - mesh.nhalo) |> abs
end

"""

Get the position of the centroid for the given _cell_ index. **NOTE:** these indices
are consistent with halo cells included. This means that if your grid has 2 halo cells,
the position of the first non-halo centroid is index at coord(mesh, 3).
The `CurvilinearGrid` only keeps track of the number of halo cells for each dimension,
whereas the grid functions have no knowledge halos. Therefore, the `coord` function
applies a shift to the index for you.
"""
@inline centroid(mesh, CI::CartesianIndex, t::Real=0) = coord(mesh, CI.I .+ 0.5, t)
@inline centroid(mesh, i::Real, t::Real=0) = coord(mesh, (i .+ 0.5,), t)

@inline function centroid(mesh, (i, j)::NTuple{2,Int}, t::Real=0)
  return coord(mesh, (i, j) .+ 0.5, t)
end

@inline function centroid(mesh, (i, j, k)::NTuple{3,Int}, t::Real=0)
  return coord(mesh, (i, j, k) .+ 0.5, t)
end

"""
Get the Jacobian matrix of the forward transformation (ξ,η,ζ) → (x,y,z).
"""
@inline jacobian_matrix(mesh, CI::CartesianIndex, t::Real=0) =
  jacobian_matrix(mesh, CI.I, t)
@inline jacobian_matrix(mesh, i::Real, t::Real=0) = jacobian_matrix(mesh, (i,), t)

"""
Get the Jacobian of the forward transformation (ξ,η,ζ) → (x,y,z).
"""
@inline jacobian(mesh, CI::CartesianIndex, t::Real=0) = jacobian(mesh, CI.I, t)
@inline jacobian(mesh, i::Real, t::Real=0) = jacobian(mesh, (i,), t)

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

function check_for_invalid_metrics(mesh::AbstractCurvilinearGrid)
  domain = mesh.iterators.cell.domain

  # cell centroid metrics
  invalid_cell_metrics = false
  for idx in domain
    if !(isfinite(mesh.cell_center_metrics.J[idx]))
      @error("Invalid jacobian @ index $(idx)")
      invalid_cell_metrics = true
    end
  end

  for (name, data) in pairs(mesh.cell_center_metrics)
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
  edge_iterators = ntuple(i -> expand_lower(domain, i, 1), length(mesh.nnodes))

  # edge metrics
  invalid_edge_metrics = false
  for (edge, edge_indices) in zip(mesh.edge_metrics, edge_iterators)
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
