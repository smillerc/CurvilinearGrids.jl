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

function update_metrics!(scheme, centroid_coordinates, cell_metrics, edge_metrics, domain)
  update_cell_center_metrics!(scheme, centroid_coordinates, cell_metrics, domain)
  update_edge_metrics!(scheme, cell_metrics, edge_metrics, domain)
  return nothing
end

function update_cell_center_metrics!(scheme, centroid_coordinates, cell_metrics, domain)
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

function update_edge_metrics!(scheme, cell_metrics, edge_metrics, domain)
  for (axis, edge) in enumerate(keys(edge_metrics)) # axis -> (1, 2, 3), edge -> (i₊½, j₊½, k₊½)
    for metric in keys(edge_metrics[edge]) # ξ, η, ζ, ξ̂, η̂, ζ̂
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
        )
      end
    end
  end

  return nothing
end

end