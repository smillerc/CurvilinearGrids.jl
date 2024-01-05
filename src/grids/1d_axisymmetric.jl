
struct CylindricalGrid1D{T1,T2,T3,T4} <: AbstractCurvilinearGrid
  r::T1
  θ::T2
  z::T3
  jacobian_matrix_func::T3
  nhalo::Int
  nnodes::Int
  limits::NamedTuple{(:ilo, :ihi),NTuple{2,Int}}
end

struct SphericalGrid1D{T1,T2,T3,T4} <: AbstractCurvilinearGrid
  r::T1
  θ::T2
  ϕ::T3
  jacobian_matrix_func::T3
  nhalo::Int
  nnodes::Int
  limits::NamedTuple{(:ilo, :ihi),NTuple{2,Int}}
end