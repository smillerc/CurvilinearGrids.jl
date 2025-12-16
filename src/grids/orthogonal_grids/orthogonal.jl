
include("cartesian_1d.jl")
include("cylindrical_1d.jl")
include("spherical_1d.jl")
include("axisymmetric_2d.jl")
include("spherical_3d.jl")
include("pad_with_halo.jl")

function _prepare_1d_coordinates(_x, nhalo, halo_coords_included)
  if !halo_coords_included
    x = pad_with_halo(_x, nhalo)
    halo_coords_included = true
  else
    x = _x
  end

  limits, iters = get_iterators(size(x), halo_coords_included, nhalo)
  celldims = size(iters.cell.full)
  nodedims = size(iters.node.full)

  return x, limits, iters, nodedims, celldims
end

function _prepare_nd_coordinates(_coords::NTuple{2}, nhalo, halo_coords_included)
  if !halo_coords_included
    coords = map(c -> pad_with_halo(c, nhalo), _coords)
    halo_coords_included = true
  else
    coords = _coords
  end

  nodedims = map(length, coords) |> Tuple
  limits, iters = get_iterators(nodedims, halo_coords_included, nhalo)
  celldims = size(iters.cell.full)

  return coords, limits, iters, nodedims, celldims
end

function _populate_1d_nodes!(storage, coords, iters, halo_coords_included)
  if halo_coords_included
    copy!(storage, coords)
  else
    @views copy!(storage[iters.node.domain], coords)
  end
  return nothing
end
