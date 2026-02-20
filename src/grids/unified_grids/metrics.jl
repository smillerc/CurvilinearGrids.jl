#
# Unified metric AoS type
#

"""
    Metric{N,T,M}

Unified AoS metric payload storing:
- `jacobian_matrix`: metric tensor (forward, inverse, or normalized inverse, depending on cache)
- `J`: determinant of `jacobian_matrix` for forward tensors, or the associated Jacobian determinant for inverse forms.
"""
struct Metric{N,T,M<:StaticMatrix{N,N,T}} 
  jacobian_matrix::M
  J::T
end

@inline function Metric(jacobian_matrix::StaticMatrix{N,N,T}) where {N,T}
  Metric{N,T,typeof(jacobian_matrix)}(jacobian_matrix, det(jacobian_matrix))
end

@inline function _metric_from_jacobian(
  jacobian_matrix::StaticMatrix{N,N,T}, J::T
) where {N,T}
  Metric{N,T,typeof(jacobian_matrix)}(jacobian_matrix, J)
end

@inline Base.getindex(metric::Metric, i::Int, j::Int) = metric.jacobian_matrix[i, j]

@inline function _metric_eltype(::Val{N}, ::Type{T}) where {N,T}
  Metric{N,T,SMatrix{N,N,T,N * N}}
end

Base.zero(::Type{Metric{N,T,M}}) where {N,T,M<:StaticMatrix{N,N,T}} = Metric(zero(M), zero(T))
Base.zero(metric::Metric{N,T,M}) where {N,T,M<:StaticMatrix{N,N,T}} = zero(Metric{N,T,M})

@inline function _as_smatrix(::Val{1}, x)
  if x isa Number
    return @SMatrix [x]
  end
  return SMatrix{1,1}(Tuple(x))
end

@inline function _as_smatrix(::Val{2}, x)
  if x isa StaticMatrix
    return SMatrix{2,2}(Tuple(x))
  end
  return SMatrix{2,2}(Tuple(x))
end

@inline function _as_smatrix(::Val{3}, x)
  if x isa StaticMatrix
    return SMatrix{3,3}(Tuple(x))
  end
  return SMatrix{3,3}(Tuple(x))
end
