
abstract type AbstractOrthogonalGrid <: AbstractCurvilinearGrid end

"""
    OrthogonalGrid{N,T,CS,...}

Unified concrete storage for orthogonal finite-volume grids across all supported
coordinate systems and dimensions.
"""
struct OrthogonalGrid{N,T,CS<:CoordinateSystemTrait,NC,CC,CV,I,DL,FA} <:
       AbstractOrthogonalGrid
  node_coordinates::NC
  centroid_coordinates::CC
  cell_volumes::CV
  iterators::I
  domain_limits::DL
  face_areas::FA
  nhalo::Int
end

# Dispatch migration guide (legacy -> unified orthogonal dispatch):
#   Old: f(g::CylindricalOrthogonalGrid1D)
#   New: f(g::OrthogonalGrid{1,<:Any,CylindricalCS})
#
#   Old: f(g::SphericalOrthogonalGrid1D)
#   New: f(g::OrthogonalGrid{1,<:Any,SphericalCS})
#
#   Old: f(g::AxisymmetricOrthogonalGrid2D)
#   New: f(g::OrthogonalGrid{2,<:Any,AxisymmetricCS{:x}})
#   New: f(g::OrthogonalGrid{2,<:Any,AxisymmetricCS{:y}})
#
#   Old: f(g::SphericalGrid3D)  # orthogonal spherical constructor return type
#   New: f(g::OrthogonalGrid{3,<:Any,SphericalCS})
#
# Prefer dispatching on `AbstractOrthogonalGrid` when the implementation is
# coordinate-system agnostic.

include("cartesian_1d.jl")
include("cartesian_2d.jl")
include("cartesian_3d.jl")
include("cylindrical_1d.jl")
include("spherical_1d.jl")
include("axisymmetric_2d.jl")
include("spherical_2d.jl")
include("spherical_3d.jl")
include("pad_with_halo.jl")

function _prepare_1d_coordinates(
  _x, nhalo, halo_coords_included; padder=pad_with_halo
)
  if !halo_coords_included
    x = padder(_x, nhalo)
  else
    x = _x
  end

  limits, iters = get_iterators(size(x), true, nhalo)
  celldims = size(iters.cell.full)
  nodedims = size(iters.node.full)

  return x, limits, iters, nodedims, celldims
end

function _prepare_nd_coordinates(
  _coords::NTuple{N}, nhalo, halo_coords_included;
  padders=ntuple(_ -> pad_with_halo, Val(N)),
) where {N}
  if !halo_coords_included
    coords = ntuple(d -> padders[d](_coords[d], nhalo), Val(N))
  else
    coords = _coords
  end

  nodedims = map(length, coords) |> Tuple
  limits, iters = get_iterators(nodedims, true, nhalo)
  celldims = size(iters.cell.full)

  return coords, limits, iters, nodedims, celldims
end

function _populate_1d_nodes!(storage, coords, iters)
  if length(coords) == length(storage)
    copy!(storage, coords)
  else
    @views storage[iters.node.domain.indices[1]] .= coords
  end
  return nothing
end

@inline function _spherical_polar_centroid(θ₀::T, θ₁::T) where {T}
  μ = clamp((cos(θ₀) + cos(θ₁)) / 2, -one(T), one(T))
  canonical = acos(μ)
  midpoint = (θ₀ + θ₁) / 2
  period = T(2π)
  positive = canonical + period * round((midpoint - canonical) / period)
  negative = -canonical + period * round((midpoint + canonical) / period)
  return abs(positive - midpoint) ≤ abs(negative - midpoint) ? positive : negative
end
