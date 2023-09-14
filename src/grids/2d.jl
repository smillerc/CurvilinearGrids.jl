
struct CurvilinearMesh2D{T1,T2,T3} <: AbstractCurvilinearMesh
  x::T1
  y::T2
  metrics::T3
  nhalo::Int
  nnodes::NTuple{2,Int}
  limits::NamedTuple{(:ilo, :ihi, :jlo, :jhi),NTuple{4,Int}}
end

function CurvilinearMesh2D(x::Function, y::Function, (n_ξ, n_η), nhalo)
  metrics = _get_metrics(x, y)
  nnodes = (n_ξ, n_η)
  ni_cells = n_ξ - 1
  nj_cells = n_η - 1
  lo = nhalo + 1
  limits = (ilo=lo, ihi=ni_cells - nhalo, jlo=lo, jhi=nj_cells - nhalo)
  return CurvilinearMesh2D(x, y, metrics, nhalo, nnodes, limits)
end

xy(m::CurvilinearMesh2D, i, j) = @SVector [m.x(i, j), m.y(i, j)]
function centroid_xy(m::CurvilinearMesh2D, i, j)
  @SVector [m.x(i + 0.5, j + 0.5), m.y(i + 0.5, j + 0.5)]
end

function coords(m::CurvilinearMesh2D)
  xy = zeros(2, m.nnodes...)
  for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i, j)
      xy[2, i, j] = m.y(i, j)
    end
  end

  return xy
end

function centroids(m::CurvilinearMesh2D)
  xy = zeros(2, (m.nnodes .- 1)...)
  for j in axes(xy, 3)
    for i in axes(xy, 2)
      xy[1, i, j] = m.x(i + 0.5, j + 0.5)
      xy[2, i, j] = m.y(i + 0.5, j + 0.5)
    end
  end

  return xy
end

function metrics(m::CurvilinearMesh2D, i, j)
  ξ̂x = m.metrics.ξ̂x
  ξ̂y = m.metrics.ξ̂y
  η̂x = m.metrics.η̂x
  η̂y = m.metrics.η̂y

  return (ξ̂x(i, j), ξ̂y(i, j), η̂x(i, j), η̂y(i, j))
end

function _get_metrics(x, y)
  ∂x∂ξ(ξ, η) = derivative(ξ -> x(ξ, η), η)
  ∂x∂η(ξ, η) = derivative(η -> x(ξ, η), ξ)
  ∂y∂ξ(ξ, η) = derivative(ξ -> y(ξ, η), η)
  ∂y∂η(ξ, η) = derivative(η -> y(ξ, η), ξ)

  jacobi_matrix(ξ, η) = @SMatrix [
    ∂x∂ξ(ξ, η) ∂x∂η(ξ, η)
    ∂y∂ξ(ξ, η) ∂y∂η(ξ, η)
  ]
  J(ξ, η) = det(jacobi_matrix(ξ, η))

  inv_jacobi_matrix(ξ, η) = inv(jacobi_matrix(ξ, η))
  J⁻¹(ξ, η) = det(inv_jacobi_matrix(ξ, η))

  ξ̂x(ξ, η) = J(ξ, η) * ∂y∂η(ξ, η)
  ξ̂y(ξ, η) = -J(ξ, η) * ∂x∂η(ξ, η)
  η̂x(ξ, η) = -J(ξ, η) * ∂y∂ξ(ξ, η)
  η̂y(ξ, η) = J(ξ, η) * ∂x∂ξ(ξ, η)

  metric_tuple = (
    ∂x∂ξ=∂x∂ξ,
    ∂x∂η=∂x∂η,
    ∂y∂ξ=∂y∂ξ,
    ∂y∂η=∂y∂η,
    J=J,
    J⁻¹=J⁻¹,
    ξ̂x=ξ̂x,
    ξ̂y=ξ̂y,
    η̂x=η̂x,
    η̂y=η̂y,
  )

  return metric_tuple
end
