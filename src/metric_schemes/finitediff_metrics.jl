module FiniteDifferenceMetrics

using FiniteDifferences
using StaticArrays
using Polyester
using LinearAlgebra

function _compute_fd_metrics(coordinates, order)
  # xyn = coords(mesh)
  xn = coordinates[1]
  yn = coordinates[2]
  ∂x_∂ξ, ∂x_∂η = _fd_metrics(xn, order)
  ∂y_∂ξ, ∂y_∂η = _fd_metrics(yn, order)

  ∂ξ_∂x = similar(xn)
  ∂η_∂x = similar(xn)
  ∂ξ_∂y = similar(xn)
  ∂η_∂y = similar(xn)

  jacobian = similar(xn)

  @batch for idx in eachindex(jacobian)
    jacobian_matrix = @SMatrix [
      ∂x_∂ξ[idx] ∂y_∂ξ[idx]
      ∂x_∂η[idx] ∂y_∂η[idx]
    ]

    inv_jacobian_matrix = inv(jacobian_matrix)

    jacobian[idx] = det(jacobian_matrix)

    ∂ξ_∂x[idx] = inv_jacobian_matrix[1, 1]
    ∂ξ_∂y[idx] = inv_jacobian_matrix[1, 2]
    ∂η_∂x[idx] = inv_jacobian_matrix[2, 1]
    ∂η_∂y[idx] = inv_jacobian_matrix[2, 2]
  end

  return (; jacobian, ∂x_∂ξ, ∂x_∂η, ∂y_∂ξ, ∂y_∂η, ∂ξ_∂x, ∂η_∂x, ∂ξ_∂y, ∂η_∂y)
end

@inline function _get_fd_jacobian_matrix(mesh, (i,)::NTuple{1,Real})
  jacobian_matrix = @SMatrix [mesh.cache.∂x_∂ξ[i]]

  return jacobian_matrix
end

@inline function _get_fd_jacobian_matrix(mesh, (i, j)::NTuple{2,Real})
  ∂x_∂ξ = mesh.cache.∂x_∂ξ
  ∂y_∂ξ = mesh.cache.∂y_∂ξ
  ∂x_∂η = mesh.cache.∂x_∂η
  ∂y_∂η = mesh.cache.∂y_∂η

  jacobian_matrix = @SMatrix [
    ∂x_∂ξ[i, j] ∂y_∂ξ[i, j]
    ∂x_∂η[i, j] ∂y_∂η[i, j]
  ]

  return jacobian_matrix
end

@inline function _get_fd_jacobian_matrix(mesh, (i, j, k)::NTuple{3,Real})
  ∂x_∂ξ = mesh.cache.∂x_∂ξ
  ∂y_∂ξ = mesh.cache.∂y_∂ξ
  ∂z_∂ξ = mesh.cache.∂z_∂ξ
  ∂x_∂η = mesh.cache.∂x_∂η
  ∂y_∂η = mesh.cache.∂y_∂η
  ∂z_∂η = mesh.cache.∂z_∂η
  ∂x_∂ζ = mesh.cache.∂x_∂ζ
  ∂y_∂ζ = mesh.cache.∂y_∂ζ
  ∂z_∂ζ = mesh.cache.∂z_∂ζ

  jacobian_matrix = @SMatrix [
    ∂x_∂ξ[i, j, k] ∂y_∂ξ[i, j, k] ∂z_∂ξ[i, j, k]
    ∂x_∂η[i, j, k] ∂y_∂η[i, j, k] ∂z_∂η[i, j, k]
    ∂x_∂ζ[i, j, k] ∂y_∂ζ[i, j, k] ∂z_∂ζ[i, j, k]
  ]

  return jacobian_matrix
end

# 1D metrics
function _fd_metrics(ϕ::AbstractArray{T,1}, order) where {T}
  c_fdm = central_fdm(order, 1)

  fullCI = CartesianIndices(ϕ)
  nhalo = (length(c_fdm.coefs)) ÷ 2
  offsets = c_fdm.grid

  ∂ϕ_∂ξ = similar(ϕ)

  iCI = @view fullCI[(begin + nhalo):(end - nhalo)]

  @inbounds for idx in iCI
    i = idx.I
    i_offset = i .+ offsets

    ∂ϕ_∂ξ[idx] = sum(view(ϕ, i_offset) .* c_fdm.coefs)
  end

  # ---------------------------------------------------
  # edges terms
  # ---------------------------------------------------
  forward_edge_offsets, forward_edge_fdms = get_edge_fdms(nhalo, order, true)
  backward_edge_offsets, backward_edge_fdms = get_edge_fdms(nhalo, order, false)

  for (offset, fdm) in zip(forward_edge_offsets, forward_edge_fdms)
    edgeCI = @view fullCI[begin + nhalo - offset]

    @inbounds for idx in edgeCI
      i = idx.I
      i_offset = i .+ fdm.grid

      ∂ϕ_∂ξ[idx] = sum(view(ϕ, i_offset) .* fdm.coefs)
    end
  end

  for (offset, fdm) in zip(backward_edge_offsets, backward_edge_fdms)
    edgeCI = @view fullCI[(end - nhalo + offset)]

    @inbounds for idx in edgeCI
      i = idx.I
      i_offset = i .+ fdm.grid

      ∂ϕ_∂ξ[idx] = sum(view(ϕ, i_offset) .* fdm.coefs)
    end
  end

  return ∂ϕ_∂ξ
end

# 2D metrics
function _fd_metrics(ϕ::AbstractArray{T,2}, order) where {T}
  c_fdm = central_fdm(order, 1)

  fullCI = CartesianIndices(ϕ)
  nhalo = (length(c_fdm.coefs)) ÷ 2
  offsets = c_fdm.grid

  ∂ϕ_∂ξ = similar(ϕ)
  ∂ϕ_∂η = similar(ϕ)

  fill!(∂ϕ_∂ξ, 0)
  fill!(∂ϕ_∂η, 0)

  ϵ = 20eps(T)

  iCI = @view fullCI[(begin + nhalo):(end - nhalo), :]
  jCI = @view fullCI[:, (begin + nhalo):(end - nhalo)]

  @batch for idx in iCI
    i, j = idx.I
    i_offset = i .+ offsets

    _∂ϕ_∂ξ = sum(view(ϕ, i_offset, j) .* c_fdm.coefs)
    ∂ϕ_∂ξ[idx] = _∂ϕ_∂ξ * (abs(_∂ϕ_∂ξ) >= ϵ)
  end

  @batch for idx in jCI
    i, j = idx.I
    j_offset = j .+ offsets

    _∂ϕ_∂η = sum(view(ϕ, i, j_offset) .* c_fdm.coefs)
    ∂ϕ_∂η[idx] = _∂ϕ_∂η * (abs(_∂ϕ_∂η) >= ϵ)
  end

  # ---------------------------------------------------
  # i direction edges for ∂ξ terms
  # ---------------------------------------------------
  forward_edge_offsets, forward_edge_fdms = get_edge_fdms(nhalo, order, true)
  backward_edge_offsets, backward_edge_fdms = get_edge_fdms(nhalo, order, false)

  for (offset, fdm) in zip(forward_edge_offsets, forward_edge_fdms)
    edgeCI = @view fullCI[begin + nhalo - offset, :]

    @batch for idx in edgeCI
      i, j = idx.I
      i_offset = i .+ fdm.grid

      _∂ϕ_∂ξ = sum(view(ϕ, i_offset, j) .* fdm.coefs)
      ∂ϕ_∂ξ[idx] = _∂ϕ_∂ξ * (abs(_∂ϕ_∂ξ) >= ϵ)
    end
  end

  for (offset, fdm) in zip(backward_edge_offsets, backward_edge_fdms)
    edgeCI = @view fullCI[(end - nhalo + offset), :]

    @batch for idx in edgeCI
      i, j = idx.I
      i_offset = i .+ fdm.grid

      _∂ϕ_∂ξ = sum(view(ϕ, i_offset, j) .* fdm.coefs)
      ∂ϕ_∂ξ[idx] = _∂ϕ_∂ξ * (abs(_∂ϕ_∂ξ) >= ϵ)
    end
  end

  # ---------------------------------------------------
  # j direction edges for ∂η terms
  # ---------------------------------------------------
  for (offset, fdm) in zip(forward_edge_offsets, forward_edge_fdms)
    edgeCI = @view fullCI[:, begin + nhalo - offset]

    @batch for idx in edgeCI
      i, j = idx.I
      j_offset = j .+ fdm.grid

      _∂ϕ_∂η = sum(view(ϕ, i, j_offset) .* fdm.coefs)
      ∂ϕ_∂η[idx] = _∂ϕ_∂η * (abs(_∂ϕ_∂η) >= ϵ)
    end
  end

  for (offset, fdm) in zip(backward_edge_offsets, backward_edge_fdms)
    edgeCI = @view fullCI[:, (end - nhalo + offset)]

    @batch for idx in edgeCI
      i, j = idx.I
      j_offset = j .+ fdm.grid

      _∂ϕ_∂η = sum(view(ϕ, i, j_offset) .* fdm.coefs)
      ∂ϕ_∂η[idx] = _∂ϕ_∂η * (abs(_∂ϕ_∂η) >= ϵ)
    end
  end

  return (∂ϕ_∂ξ, ∂ϕ_∂η)
end

# 3D metrics
function _fd_metrics(ϕ::AbstractArray{T,3}, order) where {T}
  c_fdm = central_fdm(order, 1)

  fullCI = CartesianIndices(ϕ)
  nhalo = (length(c_fdm.coefs)) ÷ 2
  offsets = c_fdm.grid

  ∂ϕ_∂ξ = similar(ϕ)
  ∂ϕ_∂η = similar(ϕ)
  ∂ϕ_∂ζ = similar(ϕ)

  iCI = @view fullCI[(begin + nhalo):(end - nhalo), :, :]
  jCI = @view fullCI[:, (begin + nhalo):(end - nhalo), :]
  kCI = @view fullCI[:, :, (begin + nhalo):(end - nhalo)]

  @inbounds for idx in iCI
    i, j, k = idx.I
    i_offset = i .+ offsets

    ∂ϕ_∂ξ[idx] = sum(view(ϕ, i_offset, j, k) .* c_fdm.coefs)
  end

  @inbounds for idx in jCI
    i, j, k = idx.I
    j_offset = j .+ offsets
    ∂ϕ_∂η[idx] = sum(view(ϕ, i, j_offset, k) .* c_fdm.coefs)
  end

  @inbounds for idx in kCI
    i, j, k = idx.I
    k_offset = k .+ offsets
    ∂ϕ_∂ζ[idx] = sum(view(ϕ, i, j, k_offset) .* c_fdm.coefs)
  end

  # ---------------------------------------------------
  # i direction edges for ∂ξ terms
  # ---------------------------------------------------
  forward_edge_offsets, forward_edge_fdms = get_edge_fdms(nhalo, order, true)
  backward_edge_offsets, backward_edge_fdms = get_edge_fdms(nhalo, order, false)

  for (offset, fdm) in zip(forward_edge_offsets, forward_edge_fdms)
    edgeCI = @view fullCI[(begin + nhalo - offset), :, :]

    @inbounds for idx in edgeCI
      i, j, k = idx.I
      i_offset = i .+ fdm.grid

      ∂ϕ_∂ξ[idx] = sum(view(ϕ, i_offset, j, k) .* fdm.coefs)
    end
  end

  for (offset, fdm) in zip(backward_edge_offsets, backward_edge_fdms)
    edgeCI = @view fullCI[(end - nhalo + offset), :, :]

    @inbounds for idx in edgeCI
      i, j, k = idx.I
      i_offset = i .+ fdm.grid

      ∂ϕ_∂ξ[idx] = sum(view(ϕ, i_offset, j, k) .* fdm.coefs)
    end
  end

  # ---------------------------------------------------
  # j direction edges for ∂η terms
  # ---------------------------------------------------
  for (offset, fdm) in zip(forward_edge_offsets, forward_edge_fdms)
    edgeCI = @view fullCI[:, (begin + nhalo - offset, :)]

    @inbounds for idx in edgeCI
      i, j, k = idx.I
      j_offset = j .+ fdm.grid

      ∂ϕ_∂η[idx] = sum(view(ϕ, i, j_offset, k) .* fdm.coefs)
    end
  end

  for (offset, fdm) in zip(backward_edge_offsets, backward_edge_fdms)
    edgeCI = @view fullCI[:, (end - nhalo + offset)]

    @inbounds for idx in edgeCI
      i, j, k = idx.I
      j_offset = j .+ fdm.grid

      ∂ϕ_∂η[idx] = sum(view(ϕ, i, j_offset, k) .* fdm.coefs)
    end
  end

  # ---------------------------------------------------
  # k direction edges for ∂ζ terms
  # ---------------------------------------------------
  for (offset, fdm) in zip(forward_edge_offsets, forward_edge_fdms)
    edgeCI = @view fullCI[:, :, (begin + nhalo - offset)]

    @inbounds for idx in edgeCI
      i, j, k = idx.I
      j_offset = j .+ fdm.grid

      ∂ϕ_∂ζ[idx] = sum(view(ϕ, i, j_offset, k) .* fdm.coefs)
    end
  end

  for (offset, fdm) in zip(backward_edge_offsets, backward_edge_fdms)
    edgeCI = @view fullCI[:, :, (end - nhalo + offset)]

    @inbounds for idx in edgeCI
      i, j, k = idx.I
      k_offset = k .+ fdm.grid

      ∂ϕ_∂ζ[idx] = sum(view(ϕ, i, j, k_offset) .* fdm.coefs)
    end
  end

  return (∂ϕ_∂ξ, ∂ϕ_∂η, ∂ϕ_∂ζ)
end

# Get forward/backward finite difference coefficients for
# the stencils that are within the halo regions along the edges
# These are hybrid forward/backward methods, e.g, using
# grid points at [-2, -1, 0, 1] rather than purely forward/backward
function get_edge_fdms(nhalo, order, forward=true, deriv=1)
  fdms = Vector{FiniteDifferences.UnadaptedFiniteDifferenceMethod{order,deriv}}(
    undef, nhalo
  )

  # the FD stencil needs to shrink appropriately 
  # due to the boundary. This shifts 
  edge_offsets = collect(1:nhalo)
  stencil_offsets = -(reverse(edge_offsets) .- 1)
  for i in edge_offsets
    # the offset for the FD stencil, e.g. -2 to +4
    offset = stencil_offsets[i]
    grid = collect(offset:(offset + order - 1))
    if !forward
      grid = -reverse(grid)
    end
    fdms[i] = FiniteDifferenceMethod(grid, 1)
  end

  return (; edge_offsets, fdms)
end

end