"""
    gcl(edge_metrics, domain::CartesianIndices)

Compute the geometric conservation law (GCL) residuals from `edge_metrics` over the provided `domain`. The return value matches the dimensionality of the domain, yielding one array per spatial dimension containing the discrete divergence of the contravariant basis vectors.
"""
@inline function _is_unified_face_metrics(em)
  em isa Tuple && !isempty(em) || return false
  first_axis = first(em)
  return first_axis isa Tuple || (first_axis isa NamedTuple && hasproperty(first_axis, :conserved))
end

@inline function _conserved_face_axis(face_axis)
  if face_axis isa NamedTuple && hasproperty(face_axis, :conserved)
    return face_axis.conserved
  end
  return face_axis[3]
end

@inline function _shift_index(idx::CartesianIndex{N}, axis::Int, δ::Int) where {N}
  CartesianIndex(ntuple(d -> d == axis ? idx.I[d] + δ : idx.I[d], N))
end

function _gcl_unified_face_metrics(em, domain::CartesianIndices{1})
  c1 = _conserved_face_axis(em[1])
  sample = c1[first(domain)]
  T = typeof(sample[1, 1])
  I₁ = similar(c1, T)
  fill!(I₁, zero(T))

  for idx in domain
    prev = _shift_index(idx, 1, -1)
    I₁[idx] = c1[idx][1, 1] - c1[prev][1, 1]
  end

  return I₁
end

function _gcl_unified_face_metrics(em, domain::CartesianIndices{2})
  c1 = _conserved_face_axis(em[1])
  c2 = _conserved_face_axis(em[2])
  sample = c1[first(domain)]
  T = typeof(sample[1, 1])
  I₁ = similar(c1, T)
  I₂ = similar(c1, T)
  fill!(I₁, zero(T))
  fill!(I₂, zero(T))

  for idx in domain
    iprev = _shift_index(idx, 1, -1)
    jprev = _shift_index(idx, 2, -1)

    I₁[idx] = (
      (c1[idx][1, 1] - c1[iprev][1, 1]) + (c2[idx][2, 1] - c2[jprev][2, 1])
    )
    I₂[idx] = (
      (c1[idx][1, 2] - c1[iprev][1, 2]) + (c2[idx][2, 2] - c2[jprev][2, 2])
    )
  end

  return I₁, I₂
end

function _gcl_unified_face_metrics(em, domain::CartesianIndices{3})
  c1 = _conserved_face_axis(em[1])
  c2 = _conserved_face_axis(em[2])
  c3 = _conserved_face_axis(em[3])
  sample = c1[first(domain)]
  T = typeof(sample[1, 1])
  I₁ = similar(c1, T)
  I₂ = similar(c1, T)
  I₃ = similar(c1, T)
  fill!(I₁, zero(T))
  fill!(I₂, zero(T))
  fill!(I₃, zero(T))

  for idx in domain
    iprev = _shift_index(idx, 1, -1)
    jprev = _shift_index(idx, 2, -1)
    kprev = _shift_index(idx, 3, -1)

    I₁[idx] = (
      (c1[idx][1, 1] - c1[iprev][1, 1]) +
      (c2[idx][2, 1] - c2[jprev][2, 1]) +
      (c3[idx][3, 1] - c3[kprev][3, 1])
    )
    I₂[idx] = (
      (c1[idx][1, 2] - c1[iprev][1, 2]) +
      (c2[idx][2, 2] - c2[jprev][2, 2]) +
      (c3[idx][3, 2] - c3[kprev][3, 2])
    )
    I₃[idx] = (
      (c1[idx][1, 3] - c1[iprev][1, 3]) +
      (c2[idx][2, 3] - c2[jprev][2, 3]) +
      (c3[idx][3, 3] - c3[kprev][3, 3])
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
