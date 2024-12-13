
using MappedArrays
using UnPack

"""
    forward_metrics!(scheme::DiscretizationScheme, metrics, x::AbstractArray{T,1}) where {T}

Compute the forward metrics ∂x/∂ξ, or `x_ξ` based on the centroid coordinates `x`.
"""
function forward_metrics!(
  scheme::DiscretizationScheme, metrics, x::AbstractArray{T,1}
) where {T}
  # axes
  ξ = 1
  domain = CartesianIndices(x)

  @views begin
    @unpack x_ξ = metrics
  end

  cell_center_derivatives!(scheme, x_ξ, x, ξ, domain)

  return nothing
end

"""
    forward_metrics!(scheme::DiscretizationScheme, metrics, x::AbstractArray{T,2}, y::AbstractArray{T,2}) where {T}

Compute the forward metrics ∂xᵢ/∂ξᵢ, i.e., `x_ξ, x_η, y_ξ, y_η`, based on the centroid coordinates `x` and `y`.
"""
function forward_metrics!(
  scheme::DiscretizationScheme,
  metrics,
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  domain,
) where {T}

  # axes
  ξ, η = (1, 2)

  if size(x) != size(y)
    error("Size mismatch for the given x, y coordinate arrays, they must all be the same")
  end

  @views begin
    xξ = metrics.x₁.ξ
    xη = metrics.x₁.η
    yξ = metrics.x₂.ξ
    yη = metrics.x₂.η
  end

  cell_center_derivatives!(scheme, xξ, x, ξ, domain)
  cell_center_derivatives!(scheme, xη, x, η, domain)
  cell_center_derivatives!(scheme, yξ, y, ξ, domain)
  cell_center_derivatives!(scheme, yη, y, η, domain)

  @kernel inbounds = true function jacobian_2d_kernel!(x_ξ, x_η, y_ξ, y_η, jacobian, domain)
    I = @index(Global, Cartesian)
    idx = domain[I]
    _jacobian_matrix = @SMatrix [
      x_ξ[idx] x_η[idx]
      y_ξ[idx] y_η[idx]
    ]

    jacobian[idx] = det(_jacobian_matrix)
  end

  jacobian_2d_kernel!(scheme.backend)(
    xξ, xη, yξ, yη, metrics.J, domain; ndrange=size(domain)
  )

  return nothing
end

"""
    forward_metrics!(scheme::DiscretizationScheme, metrics, x::AbstractArray{T,3}, y::AbstractArray{T,3}, z::AbstractArray{T,3}) where {T}

Compute the forward metrics ∂xᵢ/∂ξᵢ, i.e., `x_ξ, x_η, x_ζ, y_ξ, y_η, y_ζ, z_ξ, z_η, z_ζ`, based on the centroid coordinates `x` and `y`, and `z`.

"""
function forward_metrics!(
  scheme::DiscretizationScheme,
  metrics,
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  domain,
) where {T}

  # axes
  ξ, η, ζ = (1, 2, 3)

  if !(size(x) == size(y) == size(z))
    error(
      "Size mismatch for the given x, y, z coordinate arrays, they must all be the same"
    )
  end

  @views begin
    xξ = metrics.x₁.ξ
    xη = metrics.x₁.η
    xζ = metrics.x₁.ζ

    yξ = metrics.x₂.ξ
    yη = metrics.x₂.η
    yζ = metrics.x₂.ζ

    zξ = metrics.x₃.ξ
    zη = metrics.x₃.η
    zζ = metrics.x₃.ζ
  end

  cell_center_derivatives!(scheme, xξ, x, ξ, domain)
  cell_center_derivatives!(scheme, xη, x, η, domain)
  cell_center_derivatives!(scheme, xζ, x, ζ, domain)
  cell_center_derivatives!(scheme, yξ, y, ξ, domain)
  cell_center_derivatives!(scheme, yη, y, η, domain)
  cell_center_derivatives!(scheme, yζ, y, ζ, domain)
  cell_center_derivatives!(scheme, zξ, z, ξ, domain)
  cell_center_derivatives!(scheme, zη, z, η, domain)
  cell_center_derivatives!(scheme, zζ, z, ζ, domain)

  jacobian_3d_kernel!(scheme.backend)(
    xξ, xη, xζ, yξ, yη, yζ, zξ, zη, zζ, metrics.J, domain; ndrange=size(domain)
  )

  return nothing
end

# The compute kernel to find the jacobian in 3D
@kernel inbounds = true function jacobian_3d_kernel!(
  xξ, xη, xζ, yξ, yη, yζ, zξ, zη, zζ, J, domain
)
  I = @index(Global, Cartesian)
  idx = domain[I]
  jacobian = @SMatrix [
    xξ[idx] xη[idx] xζ[idx]
    yξ[idx] yη[idx] yζ[idx]
    zξ[idx] zη[idx] zζ[idx]
  ]

  J[idx] = det(jacobian)
end

"""
    conservative_metric!(scheme::DiscretizationScheme, ξ̂x, y_η, y_ζ, z, deriv_axes, domain)

Compute the Jacobian-normalized conservative metric. This is taken from 
Thomas P, Lombard C. AIAA J 1979;17(10):1030–7 (https://doi.org/10.2514/3.61273).
The arguments are set up to mimic the formula for the metric ξ̂x, however, it can 
be applied to any other metric, e.g., ξ̂, η̂, ζ̂ for (x,y,z). The formula is
`ξ̂x = (y_η z)_ζ − (y_ζ z)_η`, which involved inner and outer derivative terms 
(derivatives are denoted by underscores). The inner derivative terms `(y_η, y_ζ)` 
are precomputed elsewhere by the discretization scheme `scheme`.
"""
function conservative_metric!(
  scheme::DiscretizationScheme,
  ξ̂x::AbstractArray{T,3},
  y_η::AbstractArray{T,3},
  y_ζ::AbstractArray{T,3},
  z::AbstractArray{T,3},
  deriv_axes,
  domain;
  ϵ=50eps(T),
) where {T}

  #
  if isnothing(scheme.cache)
    error(
      "The DiscretizationScheme must have a cache enabled to compute conservative metrics"
    )
  end

  if !(size(ξ̂x) == size(y_η) == size(y_ζ) == size(z))
    error("Size mismatch for the given fields, they must all be the same")
  end

  ζ, η = deriv_axes

  # ξ̂x = (y_η z)_ζ − (y_ζ z)_η
  ∂ζ = mappedarray(*, y_η, z)
  ∂η = mappedarray(*, y_ζ, z)

  cell_center_derivatives!(
    scheme,
    scheme.cache.outer_deriv_1,
    scheme.cache.∂²ϕ,
    scheme.cache.∂ϕ,
    ∂ζ, # (y_η z)_ζ
    ζ, # do along the ζ dimension
    domain;
    compute_gradients=true, # compute the gradients ∂²ϕ, ∂ϕ required for (y_η z)_ζ
  ) # 1st outer deriv term

  cell_center_derivatives!(
    scheme,
    scheme.cache.outer_deriv_2,
    scheme.cache.∂²ϕ,
    scheme.cache.∂ϕ,
    ∂η, # (y_ζ z)_η
    η, # do along the η dimension
    domain;
    compute_gradients=true, # compute the gradients ∂²ϕ, ∂ϕ required for (y_ζ z)_η
  ) # 2nd outer deriv term

  @. ξ̂x = scheme.cache.outer_deriv_1 - scheme.cache.outer_deriv_2
  @. ξ̂x = ξ̂x * (abs(ξ̂x) >= ϵ)
end

"""
    conservative_metrics!(scheme::DiscretizationScheme, x::AbstractArray{T,3}, y::AbstractArray{T,3}, z::AbstractArray{T,3}, forward_metrics, inverse_metrics) where {T}

Compute all the inverse metrics ∂ξ̂ᵢ/∂xᵢ using the conservative scheme from Thomas P, Lombard C. AIAA J 1979;17(10):1030–7 (https://doi.org/10.2514/3.61273) 
"""
function conservative_metrics!(
  scheme::DiscretizationScheme,
  normalized_inverse_metrics,
  inverse_metrics,
  forward_metrics,
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  domain,
) where {T}

  # axes
  ξ, η, ζ = (1, 2, 3)

  if !(size(x) == size(y) == size(z))
    error(
      "Size mismatch for the given x, y, z coordinate arrays, they must all be the same"
    )
  end

  @views begin
    ξ̂x = normalized_inverse_metrics.ξ̂.x₁
    η̂x = normalized_inverse_metrics.η̂.x₁
    ζ̂x = normalized_inverse_metrics.ζ̂.x₁

    ξ̂y = normalized_inverse_metrics.ξ̂.x₂
    η̂y = normalized_inverse_metrics.η̂.x₂
    ζ̂y = normalized_inverse_metrics.ζ̂.x₂

    ξ̂z = normalized_inverse_metrics.ξ̂.x₃
    η̂z = normalized_inverse_metrics.η̂.x₃
    ζ̂z = normalized_inverse_metrics.ζ̂.x₃

    ξx = inverse_metrics.ξ.x₁
    ηx = inverse_metrics.η.x₁
    ζx = inverse_metrics.ζ.x₁

    ξy = inverse_metrics.ξ.x₂
    ηy = inverse_metrics.η.x₂
    ζy = inverse_metrics.ζ.x₂

    ξz = inverse_metrics.ξ.x₃
    ηz = inverse_metrics.η.x₃
    ζz = inverse_metrics.ζ.x₃

    x_ξ = forward_metrics.x₁.ξ
    x_η = forward_metrics.x₁.η
    x_ζ = forward_metrics.x₁.ζ

    y_ξ = forward_metrics.x₂.ξ
    y_η = forward_metrics.x₂.η
    y_ζ = forward_metrics.x₂.ζ

    z_ξ = forward_metrics.x₃.ξ
    z_η = forward_metrics.x₃.η
    z_ζ = forward_metrics.x₃.ζ
  end

  # ξ̂x = (y_η z)_ζ − (y_ζ z)_η
  conservative_metric!(scheme, ξ̂x, y_η, y_ζ, z, (ζ, η), domain)

  # η̂x = (y_ζ z)_ξ − (y_ξ z)_ζ
  conservative_metric!(scheme, η̂x, y_ζ, y_ξ, z, (ξ, ζ), domain)

  # ζ̂x = (y_ξ z)_η − (y_η z)_ξ
  conservative_metric!(scheme, ζ̂x, y_ξ, y_η, z, (η, ξ), domain)

  # ξ̂y = (z_η x)_ζ − (z_ζ x)_η
  conservative_metric!(scheme, ξ̂y, z_η, z_ζ, x, (ζ, η), domain)

  # η̂y = (z_ζ x)_ξ − (z_ξ x)_ζ
  conservative_metric!(scheme, η̂y, z_ζ, z_ξ, x, (ξ, ζ), domain)

  # ζ̂y = (z_ξ x)_η − (z_η x)_ξ
  conservative_metric!(scheme, ζ̂y, z_ξ, z_η, x, (η, ξ), domain)

  # ξ̂z = (x_η y)_ζ − (x_ζ y)_η
  conservative_metric!(scheme, ξ̂z, x_η, x_ζ, y, (ζ, η), domain)

  # η̂z = (x_ζ y)_ξ − (x_ξ y)_ζ
  conservative_metric!(scheme, η̂z, x_ζ, x_ξ, y, (ξ, ζ), domain)

  # ζ̂z = (x_ξ y)_η − (x_η y)_ξ
  conservative_metric!(scheme, ζ̂z, x_ξ, x_η, y, (η, ξ), domain)

  @. ξx = ξ̂x / forward_metrics.J
  @. ηx = η̂x / forward_metrics.J
  @. ζx = ζ̂x / forward_metrics.J
  @. ξy = ξ̂y / forward_metrics.J
  @. ηy = η̂y / forward_metrics.J
  @. ζy = ζ̂y / forward_metrics.J
  @. ξz = ξ̂z / forward_metrics.J
  @. ηz = η̂z / forward_metrics.J
  @. ζz = ζ̂z / forward_metrics.J

  return nothing
end

"""
    conservative_metrics!(scheme::DiscretizationScheme, normalized_inverse_metrics, inverse_metrics, forward_metrics, x::AbstractArray{T,2}, y::AbstractArray{T,2}, domain) where {T}

Compute all the inverse metrics ∂ξ̂ᵢ/∂xᵢ. This 2D version is simpler than 3D, which requires the scheme from Thomas & Lombard. 
"""
function conservative_metrics!(
  scheme::DiscretizationScheme,
  normalized_inverse_metrics,
  inverse_metrics,
  forward_metrics,
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  domain,
) where {T}
  if size(x) != size(y)
    error("Size mismatch for the given x, y coordinate arrays, they must all be the same")
  end

  @views begin
    ξ̂x = normalized_inverse_metrics.ξ̂.x₁
    η̂x = normalized_inverse_metrics.η̂.x₁
    ξ̂y = normalized_inverse_metrics.ξ̂.x₂
    η̂y = normalized_inverse_metrics.η̂.x₂

    ξx = inverse_metrics.ξ.x₁
    ηx = inverse_metrics.η.x₁
    ξy = inverse_metrics.ξ.x₂
    ηy = inverse_metrics.η.x₂

    xξ = forward_metrics.x₁.ξ
    xη = forward_metrics.x₁.η
    yξ = forward_metrics.x₂.ξ
    yη = forward_metrics.x₂.η
  end

  # This kernel finds the jacobian matrix and determinant from the 
  # forward metrics. No need to use the more complicated conservative_metric!
  # function in 2D
  inverse_metric_2d_kernel!(scheme.backend)(
    xξ, xη, yξ, yη, ξx, ξy, ηx, ηy, forward_metrics.J, domain; ndrange=size(domain)
  )

  @. ξ̂x = ξx * forward_metrics.J
  @. η̂x = ηx * forward_metrics.J
  @. ξ̂y = ξy * forward_metrics.J
  @. η̂y = ηy * forward_metrics.J

  return nothing
end

# The compute kernel to find the inverse metrics in 2D
@kernel inbounds = true function inverse_metric_2d_kernel!(
  x_ξ, x_η, y_ξ, y_η, ξ_x, ξ_y, η_x, η_y, jacobian, domain
)
  I = @index(Global, Cartesian)
  idx = domain[I]

  _jacobian_matrix = @SMatrix [
    x_ξ[idx] x_η[idx]
    y_ξ[idx] y_η[idx]
  ]

  jacobian[idx] = det(_jacobian_matrix)

  _inv_jacobian = inv(_jacobian_matrix)
  ξ_x[idx] = _inv_jacobian[1, 1]
  ξ_y[idx] = _inv_jacobian[1, 2]
  η_x[idx] = _inv_jacobian[2, 1]
  η_y[idx] = _inv_jacobian[2, 2]
end

function conservative_metrics!(
  ::DiscretizationScheme,
  inverse_metrics,
  normalized_inverse_metrics,
  forward_metrics,
  x::AbstractArray{T,1},
) where {T}
  # @. inverse_metrics.ξx = inv(forward_metrics.x_ξ)
  # @. jacobian = abs(forward_metrics.x_ξ)
  # return nothing
end

"""
    symmetric_conservative_metric!(scheme::DiscretizationScheme, ξ̂x, y_η, y_ζ, z, deriv_axes, domain)

Compute the Jacobian-normalized conservative metric. This is taken from 
T. Nonomura et al. / Computers & Fluids 107 (2015) 242–255 [http://dx.doi.org/10.1016/j.compfluid.2014.09.025]
The arguments are set up to mimic the formula for the metric ξ̂x, however, it can 
be applied to any other metric, e.g., ξ̂, η̂, ζ̂ for (x,y,z). The formula is
`ξ̂x = (1/2) * ((y_η z - z_η y)_ζ − (y_ζ z - z_ζ y)_η)`, which involved inner
and outer derivative terms (derivatives are denoted by underscores). 
The inner derivative terms `(y_η, y_ζ)` are precomputed elsewhere by the 
discretization scheme `scheme`.
"""
function symmetric_conservative_metric!(
  scheme::DiscretizationScheme,
  ξ̂x::AbstractArray{T,3},
  y_η::AbstractArray{T,3},
  z_η::AbstractArray{T,3},
  y_ζ::AbstractArray{T,3},
  z_ζ::AbstractArray{T,3},
  z::AbstractArray{T,3},
  y::AbstractArray{T,3},
  deriv_axes,
  domain,
  ϵ=50eps(T),
) where {T}

  #
  if isnothing(scheme.cache)
    error(
      "The DiscretizationScheme must have a cache enabled, i.e., DiscretizationScheme(;use_cache=true, celldims=(...)) to compute symmetric conservative metrics",
    )
  end

  ζ, η = deriv_axes

  ∂ζ = mappedarray((a, b, c, d) -> a * b - c * d, y_η, z, z_η, y)
  ∂η = mappedarray((a, b, c, d) -> a * b - c * d, y_ζ, z, z_ζ, y)

  # From Nonomura et. al
  # ξ̂x =  (1/2) * [ 
  #           (y_η z - z_η y)_ζ −   <-- this is the "∂ζ" inner term
  #           (y_ζ z - z_ζ y)_η     <-- this is the "∂η" inner term
  #      ]

  cell_center_derivatives!(
    scheme,
    scheme.cache.outer_deriv_1,
    scheme.cache.∂²ϕ,
    scheme.cache.∂ϕ,
    ∂ζ,
    ζ,
    domain;
    compute_gradients=true,
  ) # 1st outer deriv term

  cell_center_derivatives!(
    scheme,
    scheme.cache.outer_deriv_2,
    scheme.cache.∂²ϕ,
    scheme.cache.∂ϕ,
    ∂η,
    η,
    domain;
    compute_gradients=true,
  ) # 2nd outer deriv term

  @. ξ̂x = (1 / 2) * (scheme.cache.outer_deriv_1 - scheme.cache.outer_deriv_2)
  @. ξ̂x = ξ̂x * (abs(ξ̂x) >= ϵ)
end

"""
    symmetric_conservative_metrics!(scheme::DiscretizationScheme, x::AbstractArray{T,3}, y::AbstractArray{T,3}, z::AbstractArray{T,3}, forward_metrics, inverse_metrics) where {T}

Compute all the inverse metrics ∂ξ̂ᵢ/∂xᵢ using the symmetric conservative scheme from Nonomura et al. / Computers & Fluids 107 (2015) 242–255 [http://dx.doi.org/10.1016/j.compfluid.2014.09.025]
"""
function symmetric_conservative_metrics!(
  scheme::DiscretizationScheme,
  normalized_inverse_metrics,
  inverse_metrics,
  forward_metrics,
  x::AbstractArray{T,3},
  y::AbstractArray{T,3},
  z::AbstractArray{T,3},
  domain,
) where {T}
  # axes
  ξ, η, ζ = (1, 2, 3)

  if !(size(x) == size(y) == size(z))
    error(
      "Size mismatch for the given x, y, z coordinate arrays, they must all be the same"
    )
  end

  @views begin
    ξ̂x = normalized_inverse_metrics.ξ̂.x₁
    η̂x = normalized_inverse_metrics.η̂.x₁
    ζ̂x = normalized_inverse_metrics.ζ̂.x₁

    ξ̂y = normalized_inverse_metrics.ξ̂.x₂
    η̂y = normalized_inverse_metrics.η̂.x₂
    ζ̂y = normalized_inverse_metrics.ζ̂.x₂

    ξ̂z = normalized_inverse_metrics.ξ̂.x₃
    η̂z = normalized_inverse_metrics.η̂.x₃
    ζ̂z = normalized_inverse_metrics.ζ̂.x₃

    ξx = inverse_metrics.ξ.x₁
    ηx = inverse_metrics.η.x₁
    ζx = inverse_metrics.ζ.x₁

    ξy = inverse_metrics.ξ.x₂
    ηy = inverse_metrics.η.x₂
    ζy = inverse_metrics.ζ.x₂

    ξz = inverse_metrics.ξ.x₃
    ηz = inverse_metrics.η.x₃
    ζz = inverse_metrics.ζ.x₃

    x_ξ = forward_metrics.x₁.ξ
    x_η = forward_metrics.x₁.η
    x_ζ = forward_metrics.x₁.ζ

    y_ξ = forward_metrics.x₂.ξ
    y_η = forward_metrics.x₂.η
    y_ζ = forward_metrics.x₂.ζ

    z_ξ = forward_metrics.x₃.ξ
    z_η = forward_metrics.x₃.η
    z_ζ = forward_metrics.x₃.ζ
  end

  # ξ̂x = (1/2) * [(y_η z - z_η y)_ζ - (y_ζ z - z_ζ y)_η]
  # ξ̂y = (1/2) * [(z_η x - x_η z)_ζ - (z_ζ x - x_ζ z)_η]
  # ξ̂z = (1/2) * [(x_η y - y_η x)_ζ - (x_ζ y - y_ζ x)_η]
  symmetric_conservative_metric!(scheme, ξ̂x, y_η, z_η, y_ζ, z_ζ, z, y, (ζ, η), domain)
  symmetric_conservative_metric!(scheme, ξ̂y, z_η, x_η, z_ζ, x_ζ, x, z, (ζ, η), domain)
  symmetric_conservative_metric!(scheme, ξ̂z, x_η, y_η, x_ζ, y_ζ, y, x, (ζ, η), domain)

  # η̂x = (1/2) * [(y_ζ z - z_ζ y)_ξ - (y_ξ z - z_ξ y)_ζ]
  # η̂y = (1/2) * [(z_ζ x - x_ζ z)_ξ - (z_ξ x - x_ξ z)_ζ]
  # η̂z = (1/2) * [(x_ζ y - y_ζ x)_ξ - (x_ξ y - y_ξ x)_ζ]
  symmetric_conservative_metric!(scheme, η̂x, y_ζ, z_ζ, y_ξ, z_ξ, z, y, (ξ, ζ), domain)
  symmetric_conservative_metric!(scheme, η̂y, z_ζ, x_ζ, z_ξ, x_ξ, x, z, (ξ, ζ), domain)
  symmetric_conservative_metric!(scheme, η̂z, x_ζ, y_ζ, x_ξ, y_ξ, y, x, (ξ, ζ), domain)

  # ζ̂x = (1/2) * [(y_ξ z - z_ξ y)_η - (y_η z - z_η y)_ξ]
  # ζ̂y = (1/2) * [(z_ξ x - x_ξ z)_η - (z_η x - x_η z)_ξ]
  # ζ̂z = (1/2) * [(x_ξ y - y_ξ x)_η - (x_η y - y_η x)_ξ]
  symmetric_conservative_metric!(scheme, ζ̂x, y_ξ, z_ξ, y_η, z_η, z, y, (η, ξ), domain)
  symmetric_conservative_metric!(scheme, ζ̂y, z_ξ, x_ξ, z_η, x_η, x, z, (η, ξ), domain)
  symmetric_conservative_metric!(scheme, ζ̂z, x_ξ, y_ξ, x_η, y_η, y, x, (η, ξ), domain)

  @. ξx = ξ̂x / forward_metrics.J
  @. ηx = η̂x / forward_metrics.J
  @. ζx = ζ̂x / forward_metrics.J
  @. ξy = ξ̂y / forward_metrics.J
  @. ηy = η̂y / forward_metrics.J
  @. ζy = ζ̂y / forward_metrics.J
  @. ξz = ξ̂z / forward_metrics.J
  @. ηz = η̂z / forward_metrics.J
  @. ζz = ζ̂z / forward_metrics.J

  return nothing
end