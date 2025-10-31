function gcl(
  em, # edge metrics
  domain::CartesianIndices{2},
  ϵ,
)
  I₁_passes = true
  I₂_passes = true

  I₁ = -Inf
  I₂ = -Inf
  for idx in domain
    i, j = idx.I
    _I₁ = (
      (em.i₊½.ξ̂.x₁[i, j] - em.i₊½.ξ̂.x₁[i - 1, j]) +
      (em.j₊½.η̂.x₁[i, j] - em.j₊½.η̂.x₁[i, j - 1])
    )
    _I₂ = (
      (em.i₊½.ξ̂.x₂[i, j] - em.i₊½.ξ̂.x₂[i - 1, j]) +
      (em.j₊½.η̂.x₂[i, j] - em.j₊½.η̂.x₂[i, j - 1])
    )

    I₁ = max(I₁, abs(_I₁))
    I₂ = max(I₂, abs(_I₂))
  end

  I₁_passes = abs(I₁) < ϵ
  I₂_passes = abs(I₂) < ϵ

  return (I₁_passes, I₂_passes), (I₁, I₂)
end

function gcl(
  em, # edge metrics
  domain::CartesianIndices{3},
)
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