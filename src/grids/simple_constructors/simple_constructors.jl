
using KernelAbstractions
using CartesianDomains: tile

include("rectlinear.jl")
include("partitioned.jl")
include("partitioned_rectlinear.jl")
include("rtheta.jl")
include("rthetaphi.jl")

function get_subdomain_limits(cell_domain, partition_fraction, rank)
  cell_tiles = tile(cell_domain, partition_fraction; cell_based=true)
  node_tiles = tile(cell_domain, partition_fraction; cell_based=false)

  tiled_node_limits = ntuple(i -> node_tiles[i].indices, length(node_tiles))

  cell_subdomain = cell_tiles[rank] # subdomain of the cell indices
  node_subdomain = expand_upper(cell_subdomain, 1)
  on_bc = subdomain_on_bc(cell_subdomain, cell_domain)

  return on_bc, tiled_node_limits, node_subdomain, cell_subdomain
end
