using CurvilinearGrids
using CartesianDomains
using StaticArrays
using BenchmarkTools

function gcl(mesh::AbstractCurvilinearGrid2D, ϵ=5e-13)
  I₁_passes = true
  I₂_passes = true

  em = mesh.edge_metrics.inverse_normalized

  for idx in mesh.iterators.cell.domain
    i, j = idx.I
    I₁ = (
      (em.i₊½.ξ̂.x₁[i, j] - em.i₊½.ξ̂.x₁[i - 1, j]) +
      (em.j₊½.η̂.x₁[i, j] - em.j₊½.η̂.x₁[i, j - 1])
    )
    I₂ = (
      (em.i₊½.ξ̂.x₂[i, j] - em.i₊½.ξ̂.x₂[i - 1, j]) +
      (em.j₊½.η̂.x₂[i, j] - em.j₊½.η̂.x₂[i, j - 1])
    )

    I₁_passes = abs(I₁) < ϵ
    I₂_passes = abs(I₂) < ϵ
    if !(I₁_passes && I₂_passes)
      @show I₁ I₂
      break
    end
  end
  @test I₁_passes
  @test I₂_passes

  return nothing
end

function gcl(mesh::CurvilinearGrid3D, ϵ=5e-13)
  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  em = mesh.edge_metrics.inverse_normalized

  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (em.i₊½.ξ̂.x₁[i, j, k] - em.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (em.j₊½.η̂.x₁[i, j, k] - em.j₊½.η̂.x₁[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₁[i, j, k] - em.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂ = (
      (em.i₊½.ξ̂.x₂[i, j, k] - em.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (em.j₊½.η̂.x₂[i, j, k] - em.j₊½.η̂.x₂[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₂[i, j, k] - em.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃ = (
      (em.i₊½.ξ̂.x₃[i, j, k] - em.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (em.j₊½.η̂.x₃[i, j, k] - em.j₊½.η̂.x₃[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₃[i, j, k] - em.k₊½.ζ̂.x₃[i, j, k - 1])
    )

    I₁_passes = abs(I₁) < ϵ
    I₂_passes = abs(I₂) < ϵ
    I₃_passes = abs(I₃) < ϵ
    if !(I₁_passes && I₂_passes && I₃_passes)
      @show I₁ I₂ I₃
      break
    end
  end
  @test I₁_passes
  @test I₂_passes
  @test I₃_passes

  return nothing
end
