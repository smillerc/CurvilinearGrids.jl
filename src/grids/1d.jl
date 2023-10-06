struct CurvilinearGrid1D{T,T1,T2} <: AbstractCurvilinearGrid
  x::T1 # x(ξ)
  ∂x∂ξ::T2 # ∂x∂ξ(ξ)
  nhalo::Int # number of halo cells (for all dimensions)
  nnodes::NTuple{1,Int}
  limits::NamedTuple{(:ilo, :ihi),NTuple{2,Int}}
  cell_center_metrics::Vector{Metrics1D{T}}
  edge_metrics::NamedTuple{(:ξ,),NTuple{
    1, # only 1 dimension (ξ)
    Vector{Metrics1D{T}}, # an array to hold the metrics (ξx, ξt)
  }}

  J::Vector{T} # cell-centered Jacobian J
end

function CurvilinearGrid1D(x::Function, (n_ξ,), nhalo)
  ∂x∂ξ(ξ) = ForwardDiff.derivative(x, ξ)

  nnodes = (n_ξ,)
  ni_cells = n_ξ - 1
  lo = nhalo + 1
  cell_dims = (ni_cells,)

  limits = (ilo=lo, ihi=ni_cells - nhalo)

  cell_center_metrics = _make_empty_metric_array(cell_dims)

  # metric terms for the i±1/2 edges
  edge_metrics = (ξ=_make_empty_metric_array(cell_dims),)

  J = zeros(cell_dims)
  m = CurvilinearGrid1D(
    x, ∂x∂ξ, nhalo, nnodes, limits, cell_center_metrics, edge_metrics, J
  )

  # initialize the metric terms
  update(m)

  return m
end

function update(m::CurvilinearGrid1D)
  @inline for i in eachindex(m.J)
    metric = (ξ̂x=m.∂x∂ξ(i), ξt=0.0)
    edge_metric = (ξ̂x=m.∂x∂ξ(i + 0.5), ξt=0.0)
    m.cell_center_metrics[i] = metric
    m.edge_metrics.ξ[i] = edge_metric
  end

  @inbounds for i in eachindex(m.J)
    m.J[i] = m.∂x∂ξ(i)
  end

  return nothing
end

function centroids(m::CurvilinearGrid1D)
  x = zeros(m.nnodes .- 1)
  for i in axes(x, 1)
    x[i] = m.x(i + 0.5)
  end

  return x
end

function coords(m::CurvilinearGrid1D)
  x = zeros(m.nnodes)
  for i in axes(x, 1)
    x[i] = m.x(i)
  end

  return x
end
