using CurvilinearGrids
using CartesianDomains
using StaticArrays
using BenchmarkTools

# function gcl(mesh::AbstractCurvilinearGrid2D, ϵ=5e-13)
#   I₁_passes = true
#   I₂_passes = true

#   em = mesh.edge_metrics

#   for idx in mesh.iterators.cell.domain
#     i, j = idx.I
#     I₁ = (
#       (em.i₊½.ξ̂.x₁[i, j] - em.i₊½.ξ̂.x₁[i - 1, j]) +
#       (em.j₊½.η̂.x₁[i, j] - em.j₊½.η̂.x₁[i, j - 1])
#     )
#     I₂ = (
#       (em.i₊½.ξ̂.x₂[i, j] - em.i₊½.ξ̂.x₂[i - 1, j]) +
#       (em.j₊½.η̂.x₂[i, j] - em.j₊½.η̂.x₂[i, j - 1])
#     )

#     I₁_passes = abs(I₁) < ϵ
#     I₂_passes = abs(I₂) < ϵ
#     if !(I₁_passes && I₂_passes)
#       @show I₁ I₂
#       break
#     end
#   end
#   @test I₁_passes
#   @test I₂_passes

#   return nothing
# end

# function gcl(mesh::AbstractCurvilinearGrid3D, ϵ=5e-13)
#   I₁_passes = true
#   I₂_passes = true
#   I₃_passes = true

#   em = mesh.edge_metrics

#   for idx in mesh.iterators.cell.domain
#     i, j, k = idx.I
#     I₁ = (
#       (em.i₊½.ξ̂.x₁[i, j, k] - em.i₊½.ξ̂.x₁[i - 1, j, k]) +
#       (em.j₊½.η̂.x₁[i, j, k] - em.j₊½.η̂.x₁[i, j - 1, k]) +
#       (em.k₊½.ζ̂.x₁[i, j, k] - em.k₊½.ζ̂.x₁[i, j, k - 1])
#     )
#     I₂ = (
#       (em.i₊½.ξ̂.x₂[i, j, k] - em.i₊½.ξ̂.x₂[i - 1, j, k]) +
#       (em.j₊½.η̂.x₂[i, j, k] - em.j₊½.η̂.x₂[i, j - 1, k]) +
#       (em.k₊½.ζ̂.x₂[i, j, k] - em.k₊½.ζ̂.x₂[i, j, k - 1])
#     )
#     I₃ = (
#       (em.i₊½.ξ̂.x₃[i, j, k] - em.i₊½.ξ̂.x₃[i - 1, j, k]) +
#       (em.j₊½.η̂.x₃[i, j, k] - em.j₊½.η̂.x₃[i, j - 1, k]) +
#       (em.k₊½.ζ̂.x₃[i, j, k] - em.k₊½.ζ̂.x₃[i, j, k - 1])
#     )

#     I₁_passes = abs(I₁) < ϵ
#     I₂_passes = abs(I₂) < ϵ
#     I₃_passes = abs(I₃) < ϵ
#     if !(I₁_passes && I₂_passes && I₃_passes)
#       @show I₁ I₂ I₃
#       break
#     end
#   end
#   @test I₁_passes
#   @test I₂_passes
#   @test I₃_passes

#   return nothing
# end

function gcl(
  em, # edge metrics
  domain::CartesianIndices{2},
  ϵ,
)
  I₁_passes = true
  I₂_passes = true

  max_I₁ = -Inf
  max_I₂ = -Inf

  for idx in domain
    i, j = idx.I
    idx_i_prev = CartesianIndex(i - 1, j)
    idx_j_prev = CartesianIndex(i, j - 1)

    _I₁ = (
      (_hatted_metric(em, 1, 1, 1, idx) - _hatted_metric(em, 1, 1, 1, idx_i_prev)) +
      (_hatted_metric(em, 2, 2, 1, idx) - _hatted_metric(em, 2, 2, 1, idx_j_prev))
    )
    _I₂ = (
      (_hatted_metric(em, 1, 1, 2, idx) - _hatted_metric(em, 1, 1, 2, idx_i_prev)) +
      (_hatted_metric(em, 2, 2, 2, idx) - _hatted_metric(em, 2, 2, 2, idx_j_prev))
    )

    max_I₁ = max(max_I₁, abs(_I₁))
    max_I₂ = max(max_I₂, abs(_I₂))

    I₁_passes &= abs(_I₁) < ϵ
    I₂_passes &= abs(_I₂) < ϵ
  end

  return (I₁_passes, I₂_passes), (max_I₁, max_I₂)
end

function gcl(
  em, # edge metrics
  domain::CartesianIndices{3},
  ϵ,
)
  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  max_I₁ = -Inf
  max_I₂ = -Inf
  max_I₃ = -Inf

  for idx in domain
    i, j, k = idx.I
    idx_i_prev = CartesianIndex(i - 1, j, k)
    idx_j_prev = CartesianIndex(i, j - 1, k)
    idx_k_prev = CartesianIndex(i, j, k - 1)

    _I₁ = (
      (_hatted_metric(em, 1, 1, 1, idx) - _hatted_metric(em, 1, 1, 1, idx_i_prev)) +
      (_hatted_metric(em, 2, 2, 1, idx) - _hatted_metric(em, 2, 2, 1, idx_j_prev)) +
      (_hatted_metric(em, 3, 3, 1, idx) - _hatted_metric(em, 3, 3, 1, idx_k_prev))
    )
    _I₂ = (
      (_hatted_metric(em, 1, 1, 2, idx) - _hatted_metric(em, 1, 1, 2, idx_i_prev)) +
      (_hatted_metric(em, 2, 2, 2, idx) - _hatted_metric(em, 2, 2, 2, idx_j_prev)) +
      (_hatted_metric(em, 3, 3, 2, idx) - _hatted_metric(em, 3, 3, 2, idx_k_prev))
    )
    _I₃ = (
      (_hatted_metric(em, 1, 1, 3, idx) - _hatted_metric(em, 1, 1, 3, idx_i_prev)) +
      (_hatted_metric(em, 2, 2, 3, idx) - _hatted_metric(em, 2, 2, 3, idx_j_prev)) +
      (_hatted_metric(em, 3, 3, 3, idx) - _hatted_metric(em, 3, 3, 3, idx_k_prev))
    )

    max_I₁ = max(max_I₁, abs(_I₁))
    max_I₂ = max(max_I₂, abs(_I₂))
    max_I₃ = max(max_I₃, abs(_I₃))

    I₁_passes &= abs(_I₁) < ϵ
    I₂_passes &= abs(_I₂) < ϵ
    I₃_passes &= abs(_I₃) < ϵ
  end

  return (I₁_passes, I₂_passes, I₃_passes), (max_I₁, max_I₂, max_I₃)
end

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

@inline function _legacy_hatted(edge, row::Int)
  row == 1 && return edge.ξ̂
  row == 2 && return edge.η̂
  return edge.ζ̂
end

@inline function _legacy_component(hatted, comp::Int, idx::CartesianIndex)
  comp == 1 && return hatted.x₁[idx]
  comp == 2 && return hatted.x₂[idx]
  return hatted.x₃[idx]
end

@inline function _hatted_metric(em, axis::Int, row::Int, comp::Int, idx::CartesianIndex)
  if _is_unified_face_metrics(em)
    return _conserved_face_axis(em[axis])[idx][row, comp]
  end

  edge = axis == 1 ? em.i₊½ : axis == 2 ? em.j₊½ : em.k₊½
  hatted = _legacy_hatted(edge, row)
  return _legacy_component(hatted, comp, idx)
end
