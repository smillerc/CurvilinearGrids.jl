
function toedge!(
  xᵢ₊½::AbstractArray{T,N},
  ∂²x::AbstractArray{T,N},
  ∂x::AbstractArray{T,N},
  x::AbstractArray{T,N},
  domain,
  axis::Int,
  ϵ=10eps(T),
) where {T,N}

  # Compute the derivatives first
  # fill!(∂x, 0)
  # fill!(∂²x, 0)
  ∂!(∂x, x, domain, axis)
  ∂²!(∂²x, ∂x, x, domain, axis)

  # coefficients; T() is used to convert to the appropriate datatype,
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!

  a = T(1 / 2)
  b = T(1 / 12)

  # do edges
  toedge_lobc!(xᵢ₊½, ∂²x, ∂x, x, domain, axis)
  toedge_hibc!(xᵢ₊½, ∂²x, ∂x, x, domain, axis)

  inner_domain = expand(domain, axis, -1)
  for i in inner_domain
    ᵢ₊₁ = up(i, axis, 1) # i.e. [i+1, j, k], but could be i+1, j+1, or k+1
    xᴸᵢ₊½ = x[i] + a * ∂x[i] + b * ∂²x[i]
    xᴿᵢ₊½ = x[ᵢ₊₁] - a * ∂x[ᵢ₊₁] + b * ∂²x[ᵢ₊₁]

    # xᵢ₊½[i] = a * (xᴿᵢ₊½ + xᴸᵢ₊½)
    _xᵢ₊½ = xᴿᵢ₊½ + xᴸᵢ₊½
    xᵢ₊½[i] = a * (_xᵢ₊½ * (abs(_xᵢ₊½) >= ϵ))
  end

  return nothing
end

function toedge_lobc!(
  xᵢ₊½::AbstractArray{T,N},
  ∂²x::AbstractArray{T,N},
  ∂x::AbstractArray{T,N},
  x::AbstractArray{T,N},
  domain,
  axis::Int,
) where {T,N}

  # coefficients; T() is used to convert to the appropriate datatype,
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!
  a = T(1 / 2)
  b = T(1 / 12)

  b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis

  # Only the i+1/2 is stored for each cell, so at the lowest boundary
  # we have to updated the ilo-1 cell
  for i in b1
    ᵢ₋₁ = down(i, axis, 1) # i.e. [i-1, j, k], but could be i-1, j-1, or k-1
    ᵢ₊₁ = up(i, axis, 1) # i.e. [i+1, j, k], but could be i+1, j+1, or k+1

    xᴸᵢ₊½ = x[i] + a * ∂x[i] + b * ∂²x[i]
    xᴿᵢ₊½ = x[ᵢ₊₁] - a * ∂x[ᵢ₊₁] + b * ∂²x[ᵢ₊₁]

    xᵢ₊½[ᵢ₋₁] = x[i] - a * ∂x[i] + b * ∂²x[i]

    xᵢ₊½[i] = 0.5(xᴸᵢ₊½ + xᴿᵢ₊½)
  end
end

function toedge_hibc!(
  xᵢ₊½::AbstractArray{T,N},
  ∂²x::AbstractArray{T,N},
  ∂x::AbstractArray{T,N},
  x::AbstractArray{T,N},
  domain,
  axis::Int,
) where {T,N}

  # coefficients; T() is used to convert to the appropriate datatype,
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!
  a = T(1 / 2)
  b = T(1 / 12)

  b1 = upper_boundary_indices(domain, axis, 0)  # last index on given boundary axis

  for i in b1
    xᵢ₊½[i] = x[i] + a * ∂x[i] + b * ∂²x[i]
    # @show xᵢ₊½[i], x[i]
  end

  return nothing
end
