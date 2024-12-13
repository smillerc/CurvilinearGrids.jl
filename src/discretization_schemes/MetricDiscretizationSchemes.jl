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
  update_edge_metrics!(scheme, cell_metrics, edge_metrics, domain)
  return nothing
end

function update_cell_center_metrics!(
  scheme, centroid_coordinates, cell_metrics, domain; use_symmetric_metric_scheme=false
)
  forward_metrics!(
    scheme, cell_metrics.forward, StructArrays.components(centroid_coordinates)..., domain
  )

  if use_symmetric_metric_scheme
    symmetric_conservative_metrics!(
      scheme,
      cell_metrics.inverse_normalized, # ∂ξ̂ᵢ/∂xᵢ
      cell_metrics.inverse, # ∂ξᵢ/∂xᵢ
      cell_metrics.forward, # ∂xᵢ/∂ξᵢ
      StructArrays.components(centroid_coordinates)...,
    )
  else
    conservative_metrics!(
      scheme,
      cell_metrics.inverse_normalized, # ∂ξ̂ᵢ/∂xᵢ
      cell_metrics.inverse, # ∂ξᵢ/∂xᵢ
      cell_metrics.forward, # ∂xᵢ/∂ξᵢ
      StructArrays.components(centroid_coordinates)...,
    )
  end

  return nothing
end

function update_edge_metrics!(scheme, cell_metrics, edge_metrics, domain)
  for (axis, edge) in enumerate(keys(edge_metrics.inverse)) # axis -> (1, 2, 3), edge -> (i₊½, j₊½, k₊½)
    for (hat_metric, metric) in zip(
      keys(edge_metrics.inverse_normalized[edge]), # ξ̂, η̂, ζ̂
      keys(edge_metrics.inverse[edge]),            # ξ, η, ζ
    )
      for (cell_component, edge_component) in zip(
        StructArrays.components(cell_metrics.inverse_normalized[hat_metric]),       # ∂ξ̂ᵢ/∂xᵢ
        StructArrays.components(edge_metrics.inverse_normalized[edge][hat_metric]), # ∂ξ̂/∂x|ᵢ₊½
      ) # x₁, x₂, x₃

        # inverse_normalized is ξ̂ (w/ the hat)
        interpolate_to_edge!(
          scheme,
          edge_component, # ∂ξ̂/∂x|ᵢ₊½
          cell_component, # ∂ξ̂ᵢ/∂xᵢ
          axis,
          domain,
        )
      end

      for (cell_component, edge_component) in zip(
        StructArrays.components(cell_metrics.inverse[metric]),       # ∂ξᵢ/∂xᵢ
        StructArrays.components(edge_metrics.inverse[edge][metric]), # ∂ξ/∂x|ᵢ₊½
      ) # x₁, x₂, x₃

        # inverse is ξ (w/o the hat)
        interpolate_to_edge!(
          scheme,
          edge_component, # ∂ξ/∂x|ᵢ₊½
          cell_component, # ∂ξᵢ/∂xᵢ
          axis,
          domain,
        )
      end
    end
  end

  return nothing
end

end