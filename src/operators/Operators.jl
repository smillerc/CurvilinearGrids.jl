module Operators

using KernelAbstractions, Polyester, StaticArrays, LinearAlgebra
using CartesianDomains
using ..GridTypes

export cell_center_derivative, edge_derivative
# export cell_center_curl, edge_curl
export cell_center_divergence, edge_divergence
export cell_center_gradient, edge_gradient

include("derivative.jl")
include("divergence.jl")
include("gradient.jl")
# include("curl.jl") # not working yet

# --------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------

@inline cell_volume(mesh::SphericalGrid3D, I::CartesianIndex{3}) =
  @inbounds mesh.cell_volumes[I]

@inline function face_area_p(mesh::SphericalGrid3D, I::CartesianIndex{3}, axis::Int)
  @inbounds if axis == 1
    return mesh.face_areas.i₊½[I]
  elseif axis == 2
    return mesh.face_areas.j₊½[I]
  else
    return mesh.face_areas.k₊½[I]
  end
end

@inline function face_area_m(mesh::SphericalGrid3D, I::CartesianIndex{3}, axis::Int)
  I₋ = CartesianDomains.shift(I, axis, -1)
  @inbounds if axis == 1
    return mesh.face_areas.i₊½[I₋]  # i-1/2
  elseif axis == 2
    return mesh.face_areas.j₊½[I₋]  # j-1/2
  else
    return mesh.face_areas.k₊½[I₋]  # k-1/2
  end
end

@inline function face_val(A, I₁::CartesianIndex{3}, I₂::CartesianIndex{3})
  @inbounds (1 / 2) * (A[I₁] + A[I₂])
end

# Effective physical spacing at edge between I and I+ê_axis
@inline function Δx_eff_edge(mesh::SphericalGrid3D, I::CartesianIndex{3}, axis::Int)
  I₊ = CartesianDomains.shift(I, axis, +1)
  V_L = cell_volume(mesh, I)
  V_R = cell_volume(mesh, I₊)
  A_n = face_area_p(mesh, I, axis)
  return (V_L + V_R) / A_n
end

function edge_on_lo_boundary(mesh::SphericalGrid3D, cell_idx::CartesianIndex{3}, axis::Int)
  # since we're using cell indexing, the i+1/2 edge on the lo boundary needs the -1
  lo = first(mesh.iterators.cell.domain.indices[axis]) - 1
  on_bc = cell_idx.I[axis] == lo
  return on_bc
end

function edge_on_hi_boundary(mesh::SphericalGrid3D, cell_idx::CartesianIndex{3}, axis::Int)
  hi = mesh.iterators.cell.domain.indices[axis] |> last
  on_bc = cell_idx.I[axis] == hi
  return on_bc
end

end