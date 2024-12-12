module MetricDiscretizationSchemes

using KernelAbstractions
using StaticArrays
using StructArrays
using LinearAlgebra

include("DiscretizationSchemes.jl")
using .DiscretizationSchemes

include("meg/MonotoneExplicitGradients.jl")
using .MontoneExplicitGradientSchemes
export MontoneExplicitGradientScheme6thOrder

# include("central/CentralSchemes.jl")
# using .CentralSchemes
# export CentralScheme

include("metrics.jl") # forward_metrics!, conservative_metrics!

export DiscretizationScheme
export update_metrics!, update_cell_center_metrics!, update_edge_metrics!

# export forward_metrics!
# export conservative_metric!, conservative_metrics!
# export symmetric_conservative_metric!, symmetric_conservative_metrics!

function update_metrics!(scheme, centroid_coordinates, cell_metrics, edge_metrics, domain)
  update_cell_center_metrics!(scheme, centroid_coordinates, cell_metrics, domain)
  #   update_edge_metrics!(scheme, cell_metrics, edge_metrics)
  return nothing
end

function update_cell_center_metrics!(
  scheme, centroid_coordinates, cell_metrics, domain; use_symmetric_metric_scheme=false
)
  forward_metrics!(
    scheme, cell_metrics.forward, StructArrays.components(centroid_coordinates)..., domain
  )
  #   conservative_metrics!(
  #     scheme,
  #     cell_metrics.inverse_normalized, # ∂ξ̂ᵢ/∂xᵢ
  #     cell_metrics.inverse, # ∂ξᵢ/∂xᵢ
  #     cell_metrics.forward, # ∂xᵢ/∂ξᵢ
  #     centroid_coordinates...,
  #     use_symmetric_metric_scheme,
  #   )

  return nothing
end

function update_edge_metrics!(scheme, cell_metrics, edge_metrics)
  for (axis, edge) in enumerate(keys(edge_metrics.forward.normalized))

    # axis -> (1, 2, 3)
    # edge -> (i₊½, j₊½, k₊½)
    interpolate_to_edge!(
      scheme,
      cell_metrics.inverse_normalized, # ∂ξ̂ᵢ/∂xᵢ
      edge_metrics.inverse_normalized[edge], # ∂ξ̂/∂x|ᵢ₊½
      axis,
    )

    interpolate_to_edge!(
      scheme,
      cell_metrics.inverse, # ∂ξᵢ/∂xᵢ
      edge_metrics.inverse[edge], # ∂ξ/∂x|ᵢ₊½
      axis,
    )
  end

  return nothing
end

end