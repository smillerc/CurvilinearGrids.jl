module MetricDiscretizationSchemes

using KernelAbstractions
using StaticArrays
using StructArrays
using LinearAlgebra

using ..DiscretizationSchemes

include("forward_metrics.jl")
include("conservative_metrics.jl")

export update_metrics!, update_cell_center_metrics!, update_edge_metrics!

"""
    update_metrics!(scheme, centroid_coordinates, cell_metrics, edge_metrics, domain, backend)

Populate both cell-centered and edge metrics for the provided discretization `scheme` over `domain` using the supplied centroid coordinates. Edge metrics are derived from the freshly updated cell metrics. The computation is executed on the specified `backend`.
"""
function update_metrics!(
  scheme, centroid_coordinates, cell_metrics, edge_metrics, domain, backend
)
  update_cell_center_metrics!(scheme, centroid_coordinates, cell_metrics, domain, backend)
  update_edge_metrics!(scheme, cell_metrics, edge_metrics, domain, backend)
  return nothing
end

"""
    update_cell_center_metrics!(scheme, centroid_coordinates, cell_metrics, domain, backend)

Update the cell-centered metric tensors for the given `scheme` and `domain`. Depending on whether `use_symmetric_conservative_metric_scheme` is enabled on the scheme, either symmetric conservative or standard conservative metrics are computed from the `centroid_coordinates` on the provided `backend`.
"""
function update_cell_center_metrics!(
  scheme, centroid_coordinates, cell_metrics, domain, backend
)
  forward_metrics!(
    scheme, cell_metrics, StructArrays.components(centroid_coordinates)..., domain
  )

  if scheme.use_symmetric_conservative_metric_scheme
    symmetric_conservative_metrics!(
      scheme, cell_metrics, StructArrays.components(centroid_coordinates)..., domain
    )
  else
    conservative_metrics!(
      scheme, cell_metrics, StructArrays.components(centroid_coordinates)..., domain
    )
  end

  return nothing
end

"""
    update_edge_metrics!(scheme, cell_metrics, edge_metrics, domain, backend)

Interpolate cell-centered metrics onto edge-centered arrays for each axis in `domain`. The routine fills geometric Jacobians as well as metric derivatives for every edge, using `interpolate_to_edge!` to map values onto the supplied `edge_metrics` buffers.
"""
function update_edge_metrics!(scheme, cell_metrics, edge_metrics, domain, backend)
  for (axis, edge) in enumerate(keys(edge_metrics)) # axis -> (1, 2, 3), edge -> (i₊½, j₊½, k₊½)
    for metric in keys(edge_metrics[edge]) # J, ξ, η, ζ, ξ̂, η̂, ζ̂
      if metric === :J
        interpolate_to_edge!(
          scheme,
          edge_metrics[edge][:J], # Ji₊½, Jj₊½, Jk₊½)
          cell_metrics[:J], # J 
          axis,
          domain,
          backend,
        )
      else
        for (cell_component, edge_component) in zip(
          StructArrays.components(cell_metrics[metric]),       # ∂ξ̂ᵢ/∂xᵢ
          StructArrays.components(edge_metrics[edge][metric]), # ∂ξ̂/∂x|ᵢ₊½
        ) # x₁, x₂, x₃
          interpolate_to_edge!(
            scheme,
            edge_component, # ∂ξ̂/∂x|ᵢ₊½
            cell_component, # ∂ξ̂ᵢ/∂xᵢ
            axis,
            domain,
            backend,
          )
        end
      end
    end
  end

  return nothing
end

end
