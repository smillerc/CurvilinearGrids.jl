module MetricDiscretizationSchemes

using Polyester
using ChunkSplitters, KernelAbstractions
using OffsetArrays, StaticArrays
using .Threads, LinearAlgebra

export update_metrics!

include("../indexing_utils.jl")
using .IndexingUtils

include("meg6/MonotoneExplicit6thOrderScheme.jl")
using .MonotoneExplicit6thOrderScheme
export MonotoneExplicit6thOrderDiscretization
export update_edge_conserved_metrics!, update_cell_center_metrics!

function update_metrics!(
  scheme, centroids, cell_center_metrics, edge_metrics, domain; do_temporal=false
)
  update_cell_center_metrics!(scheme, cell_center_metrics, centroids, domain)
  update_edge_conserved_metrics!(
    scheme, edge_metrics, cell_center_metrics, centroids, domain
  )

  if do_temporal
    update_temporal_metrics!(scheme, edge_metrics, cell_center_metrics, domain)
  end

  return nothing
end

end
