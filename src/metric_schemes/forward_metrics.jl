
using MappedArrays
using UnPack
using CartesianDomains

"""
    forward_metrics!(scheme, metrics, x::AbstractArray{T,1}) where {T}

Compute the forward metrics ∂x/∂ξ, or `x_ξ` based on the centroid coordinates `x`.
"""
function forward_metrics!(
  scheme, metrics, x::AbstractArray{T,1}, inner_cell_domain
) where {T}
  # axes
  ξ = 1

  @views begin
    xξ = metrics.x₁.ξ
  end

  cell_center_derivatives!(
    scheme, xξ, x, ξ, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true
  )

  return nothing
end

"""
    forward_metrics!(scheme, metrics, x::AbstractArray{T,2}, y::AbstractArray{T,2}) where {T}

Compute the forward metrics ∂xᵢ/∂ξᵢ, i.e., `x_ξ, x_η, y_ξ, y_η`, based on the centroid coordinates `x` and `y`.
"""
function forward_metrics!(
  scheme, metrics, x::AbstractArray{T,2}, y::AbstractArray{T,2}, inner_cell_domain
) where {T}

  # axes
  ξ, η = (1, 2)

  if size(x) != size(y)
    error("Size mismatch for the given x, y coordinate arrays, they must all be the same")
  end

  xξ = metrics.x₁.ξ
  xη = metrics.x₁.η
  yξ = metrics.x₂.ξ
  yη = metrics.x₂.η

  cell_center_derivatives!(
    scheme, xξ, x, ξ, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true
  )

  cell_center_derivatives!(
    scheme, xη, x, η, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true
  )
  cell_center_derivatives!(
    scheme, yξ, y, ξ, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true
  )
  cell_center_derivatives!(
    scheme, yη, y, η, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true
  )

  @kernel inbounds = true function jacobian_2d_kernel!(x_ξ, x_η, y_ξ, y_η, jacobian, domain)
    I = @index(Global, Linear)
    idx = domain[I]
    _jacobian_matrix = @SMatrix [
      x_ξ[idx] x_η[idx]
      y_ξ[idx] y_η[idx]
    ]

    jacobian[idx] = det(_jacobian_matrix)
  end

  backend = KernelAbstractions.get_backend(x)

  jacobian_2d_kernel!(backend)(
    xξ, xη, yξ, yη, metrics.J, inner_cell_domain; ndrange=size(inner_cell_domain)
  )

  return nothing
end

"""
    forward_metrics!(scheme, metrics, x::AbstractArray{T,3}, y::AbstractArray{T,3}, z::AbstractArray{T,3}) where {T}

Compute the forward metrics ∂xᵢ/∂ξᵢ, i.e., `x_ξ, x_η, x_ζ, y_ξ, y_η, y_ζ, z_ξ, z_η, z_ζ`, based on the centroid coordinates `x` and `y`, and `z`.

"""
function forward_metrics!(
  scheme,
  metrics,
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  inner_cell_domain,
) where {T}

  # axes
  ξaxis, ηaxis, ζaxis = (1, 2, 3)

  if !(size(x) == size(y) == size(z))
    error(
      "Size mismatch for the given x, y, z coordinate arrays, they must all be the same"
    )
  end

  xξ = metrics.x₁.ξ
  xη = metrics.x₁.η
  xζ = metrics.x₁.ζ

  yξ = metrics.x₂.ξ
  yη = metrics.x₂.η
  yζ = metrics.x₂.ζ

  zξ = metrics.x₃.ξ
  zη = metrics.x₃.η
  zζ = metrics.x₃.ζ

  #! format: off
  cell_center_derivatives!(scheme, xξ, x, ξaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂x/∂ξ
  cell_center_derivatives!(scheme, xη, x, ηaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂x/∂η
  cell_center_derivatives!(scheme, xζ, x, ζaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂x/∂ζ 

  #
  cell_center_derivatives!(scheme, yξ, y, ξaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂y/∂ξ
  cell_center_derivatives!(scheme, yη, y, ηaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂y/∂η
  cell_center_derivatives!(scheme, yζ, y, ζaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂y/∂ζ

  #
  cell_center_derivatives!(scheme, zξ, z, ξaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂z/∂ξ
  cell_center_derivatives!(scheme, zη, z, ηaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂z/∂η
  cell_center_derivatives!(scheme, zζ, z, ζaxis, inner_cell_domain; compute_gradients=true, use_one_sided_on_edges=true) # ∂z/∂ζ
  #! format: on

  backend = KernelAbstractions.get_backend(x)
  jacobian_3d_kernel!(backend)(
    metrics.J,
    xξ,
    xη,
    xζ,
    yξ,
    yη,
    yζ,
    zξ,
    zη,
    zζ,
    inner_cell_domain;
    ndrange=size(inner_cell_domain),
  )

  return nothing
end

# The compute kernel to find the jacobian in 3D
@kernel inbounds = true function jacobian_3d_kernel!(
  J, xξ, xη, xζ, yξ, yη, yζ, zξ, zη, zζ, domain
)
  I = @index(Global, Linear)
  idx = domain[I]
  jacobian = @SMatrix [
    xξ[idx] xη[idx] xζ[idx]
    yξ[idx] yη[idx] yζ[idx]
    zξ[idx] zη[idx] zζ[idx]
  ]

  J[idx] = det(jacobian)
end
