function pad_with_halo(x::AbstractVector, nhalo::Integer)
  @assert length(x) ≥ 2 "Need at least 2 points to infer spacing"
  @assert nhalo ≥ 0

  n = length(x)
  dx_lo = x[2] - x[1]
  dx_hi = x[n] - x[n - 1]

  T = eltype(x)
  xh = similar(x, n + 2 * nhalo)

  # interior
  @inbounds for i in 1:n
    xh[nhalo + i] = x[i]
  end

  # lower halo
  @inbounds for k in 1:nhalo
    xh[nhalo + 1 - k] = x[1] - k * dx_lo
  end

  # upper halo
  @inbounds for k in 1:nhalo
    xh[nhalo + n + k] = x[n] + k * dx_hi
  end

  return xh
end

function _pad_spherical_radial_with_halo(r::AbstractVector, nhalo::Integer)
  rh = pad_with_halo(r, nhalo)

  if iszero(first(r))
    length(r) ≥ nhalo + 1 || throw(
      ArgumentError(
        "spherical radial coordinates beginning at zero require at least " *
        "$(nhalo + 1) nodes for $nhalo halo layers",
      ),
    )

    @inbounds for k in 1:nhalo
      rh[nhalo + 1 - k] = -r[k + 1]
    end
  end

  return rh
end
