"""
Shared continuous-metric support used by unified grids.

This file intentionally excludes legacy `ContinuousCurvilinearGrid*` types.
"""

include("metric_cache.jl")
include("edge_interpolation.jl")
include("cell_center_derivs.jl")

function get_iterators(celldims::NTuple{N,Int}, nhalo::Int, global_cell_domain) where {N}
  cellCI = CartesianIndices(celldims .+ 2nhalo)
  nodeCI = CartesianIndices(celldims .+ 1 .+ 2nhalo)

  node = (full=nodeCI, domain=expand(nodeCI, -nhalo))
  cell = (full=cellCI, domain=expand(cellCI, -nhalo))

  if isnothing(global_cell_domain)
    global_domain = (node=node, cell=cell)
  else
    global_domain = (node=expand_upper(global_cell_domain, +1), cell=global_cell_domain)
  end
  return (; node, cell, nhalo, global_domain)
end
