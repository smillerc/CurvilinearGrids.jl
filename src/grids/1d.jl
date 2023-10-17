
"""
CurvilinearGrid1D

# Fields
 - `x`: Node function, e.g. x -> f(ξ)
 - `∂x∂ξ`: Derivative of x wrt ξ; ∂x∂ξ(ξ)
 - `nhalo`: Number of halo cells for all dims
 - `nnodes`: Number of nodes/vertices
 - `limits`: Cell loop limits based on halo cells
"""
struct CurvilinearGrid1D{F1,F2} <: AbstractCurvilinearGrid
  x::F1 # x(ξ)
  ∂x∂ξ::F2 # ∂x∂ξ(ξ)
  nhalo::Int # number of halo cells (for all dimensions)
  nnodes::NTuple{1,Int}
  limits::NamedTuple{(:ilo, :ihi),NTuple{2,Int}}
end

"""
    CurvilinearGrid1D(x::Function, (n_ξ,), nhalo)

Create a `CurvilinearGrid1D` with a function `x(ξ)` and `nhalo` halo cells. `n_ξ` is the 
total number of nodes/vertices (not including halo).
"""
function CurvilinearGrid1D(x::F1, (n_ξ,), nhalo) where {F1<:Function}
  ∂x∂ξ(ξ) = ForwardDiff.derivative(x, ξ)
  F2 = typeof(∂x∂ξ)
  nnodes = (n_ξ,)
  ni_cells = n_ξ - 1

  # limits used for looping to avoid using the halo regions
  limits = (
    ilo=nhalo + 1, # starting index
    ihi=ni_cells + nhalo, # ending index
  )

  return CurvilinearGrid1D{F1,F2}(x, ∂x∂ξ, nhalo, nnodes, limits)
end

@inline function conservative_metrics(m::CurvilinearGrid1D, i)
  ξx = m.∂x∂ξ(i - m.nhalo)
  return (ξ̂x=1 / (ξx * ξx), ξt=zero(eltype(ξx)))
end

@inline function conservative_metrics(m::CurvilinearGrid1D, i, vx)
  # don't use (i - m.nhalo), since the static metrics() does it already
  static_metrics = conservative_metrics(m::CurvilinearGrid1D, i)
  return merge(static_metrics, (ξt=-(vx * ξx),))
end

@inline function metrics(m::CurvilinearGrid1D, i)
  ξx = m.∂x∂ξ(i - m.nhalo)
  return (ξx=1 / ξx, ξt=zero(eltype(ξx)))
end

@inline function metrics(m::CurvilinearGrid1D, i, vx)
  # don't use (i - m.nhalo), since the static metrics() does it already
  static_metrics = metrics(m, i)
  return merge(static_metrics, (ξt=-(vx * ξx),))
end

@inline function jacobian_matrix(m::CurvilinearGrid1D, i)
  return checkeps(SMatrix{1,1}(m.∂x∂ξ(i - m.nhalo)))
end

"""
    centroids(mesh::CurvilinearGrid1D, T=Float64) -> Vector{Real}

Return the vector of centroid points. This does _not_ include halo regions
since the geometry can be undefined.
"""
function centroids(mesh::CurvilinearGrid1D, T=Float64)
  x = zeros(T, cellsize(mesh))
  @inbounds for i in eachindex(x)
    x[i] = mesh.x(i + 0.5)
  end

  return x
end

"""
    coords(mesh::CurvilinearGrid1D, T=Float64) -> Vector{Real}

Return the vector of coordinate points. This does _not_ include halo regions
since the geometry can be undefined.
"""
function coords(mesh::CurvilinearGrid1D, T=Float64)

  # This doesn't account for halo regions in the indexing,
  # because the x(ξ) function in the mesh doesn't either.
  x = zeros(T, cellsize(mesh) .+ 1)
  @inbounds for i in eachindex(x)
    x[i] = mesh.x(i)
  end

  return x
end
