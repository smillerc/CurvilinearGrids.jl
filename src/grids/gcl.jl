"""
    gcl(edge_metrics, domain::CartesianIndices)

Compute the geometric conservation law (GCL) residuals from `edge_metrics` over the provided `domain`. The return value matches the dimensionality of the domain, yielding one array per spatial dimension containing the discrete divergence of the contravariant basis vectors.
"""
@inline function _is_unified_face_metrics(em)
  em isa Tuple && !isempty(em) || return false
  first_axis = first(em)
  return first_axis isa Tuple ||
         (first_axis isa NamedTuple && hasproperty(first_axis, :conserved))
end

@inline function _conserved_face_axis(face_axis)
  if face_axis isa NamedTuple && hasproperty(face_axis, :conserved)
    return face_axis.conserved
  end
  return face_axis[3]
end

function _gcl_unified_face_metrics(em, domain::CartesianIndices{1})
  ξ = 1
  x = 1

  cξ = _conserved_face_axis(em[ξ])
  sample = cξ[first(domain)]
  T = typeof(sample[ξ, x])
  I₁ = similar(cξ, T)
  fill!(I₁, zero(T))

  for idx in domain
    i, = idx.I
    idx_ξ_prev = CartesianIndex(i - 1)

    metric_ξ = cξ[idx]
    metric_ξ_prev = cξ[idx_ξ_prev]

    I₁[idx] = metric_ξ[ξ, x] - metric_ξ_prev[ξ, x]
  end

  return I₁
end

function _gcl_unified_face_metrics(em, domain::CartesianIndices{2})
  ξ = 1
  η = 2
  x = 1
  y = 2

  cξ = _conserved_face_axis(em[ξ])
  cη = _conserved_face_axis(em[η])
  sample = cξ[first(domain)]
  T = typeof(sample[ξ, x])
  I₁ = similar(cξ, T)
  I₂ = similar(cξ, T)
  fill!(I₁, zero(T))
  fill!(I₂, zero(T))

  for idx in domain
    i, j = idx.I
    idx_ξ_prev = CartesianIndex(i - 1, j)
    idx_η_prev = CartesianIndex(i, j - 1)

    metric_ξ = cξ[idx]
    metric_ξ_prev = cξ[idx_ξ_prev]
    metric_η = cη[idx]
    metric_η_prev = cη[idx_η_prev]

    I₁[idx] =
      (metric_ξ[ξ, x] - metric_ξ_prev[ξ, x]) + (metric_η[η, x] - metric_η_prev[η, x])
    I₂[idx] =
      (metric_ξ[ξ, y] - metric_ξ_prev[ξ, y]) + (metric_η[η, y] - metric_η_prev[η, y])
  end

  return I₁, I₂
end

function _gcl_unified_face_metrics(em, domain::CartesianIndices{3})
  ξ = 1
  η = 2
  ζ = 3
  x = 1
  y = 2
  z = 3

  cξ = _conserved_face_axis(em[ξ])
  cη = _conserved_face_axis(em[η])
  cζ = _conserved_face_axis(em[ζ])
  sample = cξ[first(domain)]
  T = typeof(sample[ξ, x])
  I₁ = similar(cξ, T)
  I₂ = similar(cξ, T)
  I₃ = similar(cξ, T)
  fill!(I₁, zero(T))
  fill!(I₂, zero(T))
  fill!(I₃, zero(T))

  for idx in domain
    i, j, k = idx.I
    idx_ξ_prev = CartesianIndex(i - 1, j, k)
    idx_η_prev = CartesianIndex(i, j - 1, k)
    idx_ζ_prev = CartesianIndex(i, j, k - 1)

    metric_ξ = cξ[idx]
    metric_ξ_prev = cξ[idx_ξ_prev]
    metric_η = cη[idx]
    metric_η_prev = cη[idx_η_prev]
    metric_ζ = cζ[idx]
    metric_ζ_prev = cζ[idx_ζ_prev]

    I₁[idx] = (
      (metric_ξ[ξ, x] - metric_ξ_prev[ξ, x]) +
      (metric_η[η, x] - metric_η_prev[η, x]) +
      (metric_ζ[ζ, x] - metric_ζ_prev[ζ, x])
    )
    I₂[idx] = (
      (metric_ξ[ξ, y] - metric_ξ_prev[ξ, y]) +
      (metric_η[η, y] - metric_η_prev[η, y]) +
      (metric_ζ[ζ, y] - metric_ζ_prev[ζ, y])
    )
    I₃[idx] = (
      (metric_ξ[ξ, z] - metric_ξ_prev[ξ, z]) +
      (metric_η[η, z] - metric_η_prev[η, z]) +
      (metric_ζ[ζ, z] - metric_ζ_prev[ζ, z])
    )
  end

  return I₁, I₂, I₃
end

function gcl(
  em, # edge metrics
  domain::CartesianIndices{1},
)
  if _is_unified_face_metrics(em)
    return _gcl_unified_face_metrics(em, domain)
  end

  I₁ = similar(em.i₊½.ξ̂.x₁)

  fill!(I₁, 0)

  for idx in domain
    i, = idx.I
    I₁[idx] = ((em.i₊½.ξ̂.x₁[i] - em.i₊½.ξ̂.x₁[i - 1]))
  end

  return I₁
end

function gcl(
  em, # edge metrics
  domain::CartesianIndices{2},
)
  if _is_unified_face_metrics(em)
    return _gcl_unified_face_metrics(em, domain)
  end

  I₁ = similar(em.i₊½.ξ̂.x₁)
  I₂ = similar(em.i₊½.ξ̂.x₁)

  fill!(I₁, 0)
  fill!(I₂, 0)

  for idx in domain
    i, j = idx.I
    I₁[idx] = (
      (em.i₊½.ξ̂.x₁[i, j] - em.i₊½.ξ̂.x₁[i - 1, j]) +
      (em.j₊½.η̂.x₁[i, j] - em.j₊½.η̂.x₁[i, j - 1])
    )
    I₂[idx] = (
      (em.i₊½.ξ̂.x₂[i, j] - em.i₊½.ξ̂.x₂[i - 1, j]) +
      (em.j₊½.η̂.x₂[i, j] - em.j₊½.η̂.x₂[i, j - 1])
    )
  end

  return I₁, I₂
end

function gcl(
  em, # edge metrics
  domain::CartesianIndices{3},
)
  if _is_unified_face_metrics(em)
    return _gcl_unified_face_metrics(em, domain)
  end

  I₁ = similar(em.i₊½.ξ̂.x₁)
  I₂ = similar(em.i₊½.ξ̂.x₁)
  I₃ = similar(em.i₊½.ξ̂.x₁)

  fill!(I₁, 0)
  fill!(I₂, 0)
  fill!(I₃, 0)

  for idx in domain
    i, j, k = idx.I
    I₁[idx] = (
      (em.i₊½.ξ̂.x₁[i, j, k] - em.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (em.j₊½.η̂.x₁[i, j, k] - em.j₊½.η̂.x₁[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₁[i, j, k] - em.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂[idx] = (
      (em.i₊½.ξ̂.x₂[i, j, k] - em.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (em.j₊½.η̂.x₂[i, j, k] - em.j₊½.η̂.x₂[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₂[i, j, k] - em.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃[idx] = (
      (em.i₊½.ξ̂.x₃[i, j, k] - em.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (em.j₊½.η̂.x₃[i, j, k] - em.j₊½.η̂.x₃[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₃[i, j, k] - em.k₊½.ζ̂.x₃[i, j, k - 1])
    )
  end

  return I₁, I₂, I₃
end
