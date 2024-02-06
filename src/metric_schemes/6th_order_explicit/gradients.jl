function ∂!(
  ∂x::AbstractArray{T,N}, x::AbstractArray{T,N}, domain, axis::Int, ϵ=2eps(T)
) where {T,N}

  # for small domains, skip, since the deriv will be zero
  if size(domain, axis) <= 2
    fill!(∂x, 0)
    return nothing
  end

  # do edges
  ∂_at_lobc!(∂x, x, domain, axis)
  ∂_at_hibc!(∂x, x, domain, axis)

  # coefficients; T() is used to convert to the appropriate datatype, 
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!
  a = T(3 / 4)
  b = T(3 / 20)
  c = T(1 / 60)

  if 2 < size(domain, axis) <= 4 # 2nd order
    inner_domain = expand(domain, axis, -1)
    @info "2nd order"
    for i in inner_domain
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)
      _∂x = (x[ᵢ₊₁] - x[ᵢ₋₁]) / 2
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end

  elseif 4 < size(domain, axis) <= 7 # 4th order
    @info "4th order"
    inner_domain = expand(domain, axis, -2)
    for i in inner_domain
      ᵢ₋₂ = down(i, axis, 2)
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₊₂ = up(i, axis, 2)

      _∂x = (-x[ᵢ₊₂] + 8x[ᵢ₊₁] - 8x[ᵢ₋₁] + x[ᵢ₋₂]) / 12
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end

  else # 6th order; size(domain, axis) > 7
    inner_domain = expand(domain, axis, -3)
    for i in inner_domain
      ᵢ₋₃ = down(i, axis, 3)
      ᵢ₋₂ = down(i, axis, 2)
      ᵢ₋₁ = down(i, axis, 1)

      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₊₂ = up(i, axis, 2)
      ᵢ₊₃ = up(i, axis, 3)

      _∂x = (a * (x[ᵢ₊₁] - x[ᵢ₋₁]) - b * (x[ᵢ₊₂] - x[ᵢ₋₂]) + c * (x[ᵢ₊₃] - x[ᵢ₋₃]))
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end
  end

  return nothing
end

function ∂²!(
  ∂²x::AbstractArray{T,N},
  ∂x::AbstractArray{T,N},
  x::AbstractArray{T,N},
  domain,
  axis::Int,
  ϵ=2eps(T),
) where {T,N}

  # for small domains, skip, since the deriv will be zero
  if size(domain, axis) <= 2
    fill!(∂²x, 0)
    return nothing
  end

  # do edges
  ∂²_at_lobc!(∂²x, ∂x, x, domain, axis)
  ∂²_at_hibc!(∂²x, ∂x, x, domain, axis)

  inner_domain = expand(domain, axis, -1)
  # @info "∂²!", inner_domain

  for i in inner_domain
    ᵢ₋₁ = down(i, axis, 1)
    ᵢ₊₁ = up(i, axis, 1)

    _∂²x = 2(x[ᵢ₊₁] - 2x[i] + x[ᵢ₋₁]) - (∂x[ᵢ₊₁] - ∂x[ᵢ₋₁]) / 2
    ∂²x[i] = _∂²x * (abs(_∂²x) >= ϵ)
  end

  return nothing
end

function ∂_at_lobc!(
  ∂x::AbstractArray{T,N}, x::AbstractArray{T,N}, domain, axis::Int, ϵ=2eps(T)
) where {T,N}

  # coefficients; T() is used to convert to the appropriate datatype, 
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!
  a = T(1 / 4)
  b = T(5 / 6)
  c = T(3 / 2)
  d = T(1 / 2)
  e = T(1 / 12)
  f = T(2 / 3)

  if size(domain, axis) > 2
    b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis
    for i in b1
      ᵢ₊₁ = up(i, axis, +1)
      ᵢ₊₂ = up(i, axis, +2)
      _∂x = -c * x[i] + 2x[ᵢ₊₁] - d * x[ᵢ₊₂]
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end
  end

  if size(domain, axis) > 4
    b2 = lower_boundary_indices(domain, axis, +1) # first index + 1 on given boundary axis
    for i in b2
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₊₂ = up(i, axis, 2)
      ᵢ₊₃ = up(i, axis, 3)

      _∂x = -a * x[ᵢ₋₁] - b * x[i] + c * x[ᵢ₊₁] - d * x[ᵢ₊₂] + e * x[ᵢ₊₃]
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end
  end

  if size(domain, axis) > 6
    b3 = lower_boundary_indices(domain, axis, +2) # first index + 2 on given boundary axis
    for i in b3
      ᵢ₋₂ = down(i, axis, 2)
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₊₂ = up(i, axis, 2)

      _∂x = f * (x[ᵢ₊₁] - x[ᵢ₋₁]) - e * (x[ᵢ₊₂] - x[ᵢ₋₂])
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end
  end

  return nothing
end

function ∂_at_hibc!(
  ∂x::AbstractArray{T,N}, x::AbstractArray{T,N}, domain, axis::Int, ϵ=2eps(T)
) where {T,N}

  # coefficients; T() is used to convert to the appropriate datatype, 
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!
  a = T(1 / 4)
  b = T(5 / 6)
  c = T(3 / 2)
  d = T(1 / 2)
  e = T(1 / 12)
  f = T(2 / 3)

  if size(domain, axis) > 2
    b1 = upper_boundary_indices(domain, axis, 0)  # first index on given boundary axis
    for i in b1
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₋₂ = down(i, axis, 2)
      _∂x = c * x[i] - 2x[ᵢ₋₁] + d * x[ᵢ₋₂]
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end
  end

  if size(domain, axis) > 4
    b2 = upper_boundary_indices(domain, axis, -1) # first index + 1 on given boundary axis
    for i in b2
      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₋₂ = down(i, axis, 2)
      ᵢ₋₃ = down(i, axis, 3)

      _∂x = a * x[ᵢ₊₁] + b * x[i] - c * x[ᵢ₋₁] + d * x[ᵢ₋₂] - e * x[ᵢ₋₃]
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end
  end

  if size(domain, axis) > 6
    b3 = upper_boundary_indices(domain, axis, -2) # first index + 2 on given boundary axis
    for i in b3
      ᵢ₋₂ = down(i, axis, 2)
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₊₂ = up(i, axis, 2)

      _∂x = f * (x[ᵢ₊₁] - x[ᵢ₋₁]) - e * (x[ᵢ₊₂] - x[ᵢ₋₂])
      ∂x[i] = _∂x * (abs(_∂x) >= ϵ)
    end
  end

  return nothing
end

function ∂²_at_lobc!(
  ∂²x::AbstractArray{T,N},
  ∂x::AbstractArray{T,N},
  x::AbstractArray{T,N},
  domain,
  axis::Int,
  ϵ=2eps(T),
) where {T,N}

  # coefficients; T() is used to convert to the appropriate datatype, 
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!
  a = T(1 / 4)
  b = T(5 / 6)
  c = T(3 / 2)
  d = T(1 / 2)
  e = T(1 / 12)
  f = T(2 / 3)

  if size(domain, axis) > 2
    b1 = lower_boundary_indices(domain, axis, 0)  # first index on given boundary axis
    for i in b1
      ᵢ₊₁ = up(i, axis, +1)
      ᵢ₊₂ = up(i, axis, +2)
      _∂²x = -c * ∂x[i] + 2∂x[ᵢ₊₁] - d * ∂x[ᵢ₊₂]
      ∂²x[i] = _∂²x * (abs(_∂²x) >= ϵ)
    end
  end

  if size(domain, axis) > 4
    b2 = lower_boundary_indices(domain, axis, +1) # first index + 1 on given boundary axis
    for i in b2
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₊₂ = up(i, axis, 2)
      ᵢ₊₃ = up(i, axis, 3)

      _∂²x = -a * ∂x[ᵢ₋₁] - b * ∂x[i] + c * ∂x[ᵢ₊₁] - d * ∂x[ᵢ₊₂] + e * ∂x[ᵢ₊₃]
      ∂²x[i] = _∂²x * (abs(_∂²x) >= ϵ)
    end
  end

  if size(domain, axis) > 6
    b3 = lower_boundary_indices(domain, axis, +2) # first index + 2 on given boundary axis
    for i in b3
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)

      _∂²x = f * (∂x[ᵢ₊₁] - ∂x[ᵢ₋₁]) - 2(x[ᵢ₊₁] - 2x[i] + x[ᵢ₋₁])
      ∂²x[i] = _∂²x * (abs(_∂²x) >= ϵ)
    end
  end
  return nothing
end

function ∂²_at_hibc!(
  ∂²x::AbstractArray{T,N},
  ∂x::AbstractArray{T,N},
  x::AbstractArray{T,N},
  domain,
  axis::Int,
  ϵ=2eps(T),
) where {T,N}

  # coefficients; T() is used to convert to the appropriate datatype, 
  # which is more impactful on GPUs, e.g. Float32 is significantly faster.
  # These are all compile-time constants, so the cost is nothing!
  a = T(1 / 4)
  b = T(5 / 6)
  c = T(3 / 2)
  d = T(1 / 2)
  e = T(1 / 12)
  f = T(2 / 3)

  if size(domain, axis) > 2
    b1 = upper_boundary_indices(domain, axis, 0)  # last index on given boundary axis
    for i in b1
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₋₂ = down(i, axis, 2)
      _∂²x = c * ∂x[i] - 2∂x[ᵢ₋₁] + d * ∂x[ᵢ₋₂]
      ∂²x[i] = _∂²x * (abs(_∂²x) >= ϵ)
    end
  end

  if size(domain, axis) > 4
    b2 = upper_boundary_indices(domain, axis, -1) # last index - 1 on given boundary axis
    for i in b2
      ᵢ₊₁ = up(i, axis, 1)
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₋₂ = down(i, axis, 2)
      ᵢ₋₃ = down(i, axis, 3)

      _∂²x = a * ∂x[ᵢ₊₁] + b * ∂x[i] - c * ∂x[ᵢ₋₁] + d * ∂x[ᵢ₋₂] - e * ∂x[ᵢ₋₃]
      ∂²x[i] = _∂²x * (abs(_∂²x) >= ϵ)
    end
  end

  if size(domain, axis) > 6
    b3 = upper_boundary_indices(domain, axis, -2) # last index - 2 on given boundary axis
    for i in b3
      ᵢ₋₁ = down(i, axis, 1)
      ᵢ₊₁ = up(i, axis, 1)

      _∂²x = f * (∂x[ᵢ₊₁] - ∂x[ᵢ₋₁]) - 2(x[ᵢ₊₁] - 2x[i] + x[ᵢ₋₁])
      ∂²x[i] = _∂²x * (abs(_∂²x) >= ϵ)
    end
  end

  return nothing
end
