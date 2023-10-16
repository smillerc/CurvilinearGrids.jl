
"""
CurvilinearGrid1D

# Fields
 - `x` : Node function, e.g. x -> f(ξ)
 - `∂x∂ξ` : Derivative of x wrt ξ; ∂x∂ξ(ξ)
 - `nhalo` : Number of halo cells for all dims
 - `nnodes` : Number of nodes/vertices
 - `limits` : Cell loop limits based on halo cells
 - `cell_center_metrics` : Grid metrics at each cell center
 - `edge_metrics` : Grid metrics at each i+1/2 edge
 - `J` : cell-centered Jacobian

"""
struct CurvilinearGrid1D{T1,T2} <: AbstractCurvilinearGrid
  x::T1 # x(ξ)
  ∂x∂ξ::T2 # ∂x∂ξ(ξ)
  nhalo::Int # number of halo cells (for all dimensions)
  nnodes::NTuple{1,Int}
  limits::NamedTuple{(:ilo, :ihi),NTuple{2,Int}}
  # cell_center_metrics::Vector{Metrics1D{T}}
  # edge_metrics::NamedTuple{(:ξ,),NTuple{
  #   1, # only 1 dimension (ξ)
  #   Vector{Metrics1D{T}}, # an array to hold the metrics (ξx, ξt)
  # }}

  # J::Vector{T} # cell-centered Jacobian J
end

"""
    CurvilinearGrid1D(x::Function, (n_ξ,), nhalo)

Create a `CurvilinearGrid1D` with a function `x(ξ)` and `nhalo` halo cells. `n_ξ` is the 
total number of nodes/vertices (not including halo).
"""
function CurvilinearGrid1D(x::Function, (n_ξ,), nhalo)
  ∂x∂ξ(ξ) = ForwardDiff.derivative(x, ξ)

  nnodes = (n_ξ,)
  ni_cells = n_ξ - 1
  lo = nhalo + 1
  # cell_dims = (ni_cells,)
  # cell_dims_whalo = cell_dims .+ 2nhalo
  limits = (ilo=lo, ihi=ni_cells + nhalo)

  # cell_center_metrics = _make_empty_metric_array(cell_dims_whalo)

  # # metric terms for the i±1/2 edges
  # edge_metrics = (ξ=_make_empty_metric_array(cell_dims_whalo),)

  # J = zeros(cell_dims_whalo)
  m = CurvilinearGrid1D(
    x,
    ∂x∂ξ,
    nhalo,
    nnodes,
    limits,
    # cell_center_metrics, edge_metrics, J
  )

  return m
end

function conservative_metrics(m::CurvilinearGrid1D, i)
  ξx = m.∂x∂ξ(i)
  return (ξ̂x=1 / (ξx * ξx), ξt=zero(eltype(ξx)))
end

function metrics(m::CurvilinearGrid1D, i)
  ξx = m.∂x∂ξ(i)
  return (ξx=1 / ξx, ξt=zero(eltype(ξx)))
end

function metrics(m::CurvilinearGrid1D, i, vx)
  static_metrics = metrics(m, i)
  return merge(static_metrics, (ξt=-(vx * ξx),))
end

function jacobian_matrix(m::CurvilinearGrid1D, i)
  return SMatrix{1,1}(m.∂x∂ξ(i))
end

function centroids(m::CurvilinearGrid1D)
  x = zeros(cellsize(m))
  ilo, ihi = m.limits
  for i in ilo:ihi
    x[i] = m.x(i + 0.5)
  end

  return x
end

function coords(mesh::CurvilinearGrid1D)
  @unpack ilo, ihi = mesh.limits
  return mesh.x.(ilo:(ihi + 1))
end
