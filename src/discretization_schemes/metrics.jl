
using MappedArrays
using UnPack

"""
    forward_metrics!(scheme::DiscretizationScheme, metrics, x::AbstractArray{T,1}) where {T}

Compute the forward metrics ∂x/∂ξ, or `x_ξ` based on the centroid coordinates `x`.
"""
function forward_metrics!(
  scheme::DiscretizationScheme, metrics, x::AbstractArray{T,1}, domain
) where {T}
  # axes
  ξ = 1

  @views begin
    xξ = metrics.x₁.ξ
  end

  cell_center_derivatives!(scheme, xξ, x, ξ, domain)

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
  # ∂ζ = mappedarray(*, y_η, z) :: This is a mapped array, which is incompatible with GPU kernels
  scheme.cache.inner_deriv_1 .= y_η .* z
  # ∂η = mappedarray(*, y_ζ, z) :: This is a mapped array, which is incompatible with GPU kernels
  scheme.cache.inner_deriv_2 .= y_ζ .* z

  cell_center_derivatives!(
    scheme,
    scheme.cache.outer_deriv_1,
    scheme.cache.∂²ϕ,
    scheme.cache.∂ϕ,
    scheme.cache.inner_deriv_1, # (y_η z)_ζ
    ζ, # do along the ζ dimension
    domain;
    compute_gradients=true, # compute the gradients ∂²ϕ, ∂ϕ required for (y_η z)_ζ
  ) # 1st outer deriv term

  cell_center_derivatives!(
    scheme,
    scheme.cache.outer_deriv_2,
    scheme.cache.∂²ϕ,
    scheme.cache.∂ϕ,
    scheme.cache.inner_deriv_2, # (y_ζ z)_η
    η, # do along the η dimension
    domain;
    compute_gradients=true, # compute the gradients ∂²ϕ, ∂ϕ required for (y_ζ z)_η
  ) # 2nd outer deriv term

  @views begin
    @. ξ̂x[domain] = scheme.cache.outer_deriv_1[domain] - scheme.cache.outer_deriv_2[domain]
    @. ξ̂x[domain] = ξ̂x[domain] * (abs(ξ̂x[domain]) >= ϵ)
  end
end

"""
    conservative_metrics!(scheme::DiscretizationScheme, metrics, x::AbstractArray{T,3}, y::AbstractArray{T,3}, z::AbstractArray{T,3}, domain) where {T}

Compute all the inverse metrics ∂ξ̂ᵢ/∂xᵢ using the conservative scheme from Thomas P, Lombard C. AIAA J 1979;17(10):1030–7 (https://doi.org/10.2514/3.61273) 
"""
function conservative_metrics!(
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
    ξ̂x = metrics.ξ̂.x₁
    η̂x = metrics.η̂.x₁
    ζ̂x = metrics.ζ̂.x₁

    ξ̂y = metrics.ξ̂.x₂
    η̂y = metrics.η̂.x₂
    ζ̂y = metrics.ζ̂.x₂

    ξ̂z = metrics.ξ̂.x₃
    η̂z = metrics.η̂.x₃
    ζ̂z = metrics.ζ̂.x₃

    ξx = metrics.ξ.x₁
    ηx = metrics.η.x₁
    ζx = metrics.ζ.x₁

    ξy = metrics.ξ.x₂
    ηy = metrics.η.x₂
    ζy = metrics.ζ.x₂

    ξz = metrics.ξ.x₃
    ηz = metrics.η.x₃
    ζz = metrics.ζ.x₃

    x_ξ = metrics.x₁.ξ
    x_η = metrics.x₁.η
    x_ζ = metrics.x₁.ζ

    y_ξ = metrics.x₂.ξ
    y_η = metrics.x₂.η
    y_ζ = metrics.x₂.ζ

    z_ξ = metrics.x₃.ξ
    z_η = metrics.x₃.η
    z_ζ = metrics.x₃.ζ
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

  @views begin
    @. ξx[domain] = ξ̂x[domain] / metrics.J[domain]
    @. ηx[domain] = η̂x[domain] / metrics.J[domain]
    @. ζx[domain] = ζ̂x[domain] / metrics.J[domain]
    @. ξy[domain] = ξ̂y[domain] / metrics.J[domain]
    @. ηy[domain] = η̂y[domain] / metrics.J[domain]
    @. ζy[domain] = ζ̂y[domain] / metrics.J[domain]
    @. ξz[domain] = ξ̂z[domain] / metrics.J[domain]
    @. ηz[domain] = η̂z[domain] / metrics.J[domain]
    @. ζz[domain] = ζ̂z[domain] / metrics.J[domain]
  end

  return nothing
end

"""
    conservative_metrics!(scheme::DiscretizationScheme, metrics, x::AbstractArray{T,2}, y::AbstractArray{T,2}, domain) where {T}

Compute all the inverse metrics ∂ξ̂ᵢ/∂xᵢ. This 2D version is simpler than 3D, which requires the scheme from Thomas & Lombard. 
"""
function conservative_metrics!(
  scheme::DiscretizationScheme,
  metrics,
  x::AbstractArray{T,2},
  y::AbstractArray{T,2},
  domain,
) where {T}
  if size(x) != size(y)
    error("Size mismatch for the given x, y coordinate arrays, they must all be the same")
  end

  ξ̂x = metrics.ξ̂.x₁
  η̂x = metrics.η̂.x₁
  ξ̂y = metrics.ξ̂.x₂
  η̂y = metrics.η̂.x₂

  ξx = metrics.ξ.x₁
  ηx = metrics.η.x₁
  ξy = metrics.ξ.x₂
  ηy = metrics.η.x₂

  xξ = metrics.x₁.ξ
  xη = metrics.x₁.η
  yξ = metrics.x₂.ξ
  yη = metrics.x₂.η

  # This kernel finds the jacobian matrix and determinant from the 
  # forward metrics. No need to use the more complicated conservative_metric!
  # function in 2D
  inverse_metric_2d_kernel!(scheme.backend)(
    xξ, xη, yξ, yη, ξx, ξy, ηx, ηy, metrics.J, domain; ndrange=size(domain)
  )

  @views begin
    @. ξ̂x[domain] = ξx[domain] * metrics.J[domain]
    @. η̂x[domain] = ηx[domain] * metrics.J[domain]
    @. ξ̂y[domain] = ξy[domain] * metrics.J[domain]
    @. η̂y[domain] = ηy[domain] * metrics.J[domain]
  end

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

# Does this require domain limiting?
function conservative_metrics!(
  ::DiscretizationScheme, metrics, x::AbstractArray{T,1}, domain
) where {T}
  @. metrics.ξ.x₁ = inv(metrics.x₁.ξ)
  @. metrics.J = abs(metrics.x₁.ξ)
  @. metrics.ξ̂.x₁ = metrics.ξ.x₁ * metrics.J

  return nothing
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
    symmetric_conservative_metrics!(scheme::DiscretizationScheme, metrics, x::AbstractArray{T,3}, y::AbstractArray{T,3}, z::AbstractArray{T,3}, domain) where {T}

Compute all the inverse metrics ∂ξ̂ᵢ/∂xᵢ using the symmetric conservative scheme from Nonomura et al. / Computers & Fluids 107 (2015) 242–255 [http://dx.doi.org/10.1016/j.compfluid.2014.09.025]
"""
function symmetric_conservative_metrics!(
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
    ξ̂x = metrics.ξ̂.x₁
    η̂x = metrics.η̂.x₁
    ζ̂x = metrics.ζ̂.x₁

    ξ̂y = metrics.ξ̂.x₂
    η̂y = metrics.η̂.x₂
    ζ̂y = metrics.ζ̂.x₂

    ξ̂z = metrics.ξ̂.x₃
    η̂z = metrics.η̂.x₃
    ζ̂z = metrics.ζ̂.x₃

    ξx = metrics.ξ.x₁
    ηx = metrics.η.x₁
    ζx = metrics.ζ.x₁

    ξy = metrics.ξ.x₂
    ηy = metrics.η.x₂
    ζy = metrics.ζ.x₂

    ξz = metrics.ξ.x₃
    ηz = metrics.η.x₃
    ζz = metrics.ζ.x₃

    x_ξ = metrics.x₁.ξ
    x_η = metrics.x₁.η
    x_ζ = metrics.x₁.ζ

    y_ξ = metrics.x₂.ξ
    y_η = metrics.x₂.η
    y_ζ = metrics.x₂.ζ

    z_ξ = metrics.x₃.ξ
    z_η = metrics.x₃.η
    z_ζ = metrics.x₃.ζ
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

  @. ξx = ξ̂x / metrics.J
  @. ηx = η̂x / metrics.J
  @. ζx = ζ̂x / metrics.J
  @. ξy = ξ̂y / metrics.J
  @. ηy = η̂y / metrics.J
  @. ζy = ζ̂y / metrics.J
  @. ξz = ξ̂z / metrics.J
  @. ηz = η̂z / metrics.J
  @. ζz = ζ̂z / metrics.J

  return nothing
end
