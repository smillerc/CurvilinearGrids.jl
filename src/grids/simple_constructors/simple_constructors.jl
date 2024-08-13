
using KernelAbstractions
using CartesianDomains: tile

include("rectlinear.jl")
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

function subdomain_on_bc(subdomain, domain::CartesianIndices{1})
  return (
    ilo=first(subdomain.indices[1]) == first(domain.indices[1]),
    ihi=last(subdomain.indices[1]) - 1 == last(domain.indices[1]),
  )
end

function subdomain_on_bc(subdomain, domain::CartesianIndices{2})
  return (
    ilo=first(subdomain.indices[1]) == first(domain.indices[1]),
    ihi=last(subdomain.indices[1]) - 1 == last(domain.indices[1]),
    jlo=first(subdomain.indices[2]) == first(domain.indices[2]),
    jhi=last(subdomain.indices[2]) - 1 == last(domain.indices[2]),
  )
end

function subdomain_on_bc(subdomain, domain::CartesianIndices{3})
  return (
    ilo=first(subdomain.indices[1]) == first(domain.indices[1]),
    ihi=last(subdomain.indices[1]) - 1 == last(domain.indices[1]),
    jlo=first(subdomain.indices[2]) == first(domain.indices[2]),
    jhi=last(subdomain.indices[2]) - 1 == last(domain.indices[2]),
    klo=first(subdomain.indices[3]) == first(domain.indices[3]),
    khi=last(subdomain.indices[3]) - 1 == last(domain.indices[3]),
  )
end