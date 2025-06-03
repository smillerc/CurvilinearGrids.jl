module SpecialArrays

export RectlinearArray

struct RectlinearArray{T,N,D} <: AbstractArray{T,N}
    data::AbstractArray{T,D}
    dims::NTuple{N,Int}
    fixed_index::Tuple{Vararg{Int}}
end

function RectlinearArray(data::AbstractArray{T}, fixed_index::Tuple) where {T}
    new_data = _drop_dims(data, fixed_index)
    return RectlinearArray{eltype(data), ndims(data), ndims(new_data)}(new_data, size(data), fixed_index)
end

Base.size(A::RectlinearArray) = A.dims

Base.getindex(A::RectlinearArray{T,N}, I::Vararg{Int, N}) where {T,N} = _skip_index(A.data, A.fixed_index, I...)

Base.setindex!(A::RectlinearArray{T,N}, v, I::Vararg{Int,N}) where {T,N} = (A.data[_drop_index(I, A.fixed_index)...] = v)

Base.BroadcastStyle(::Type{<:RectlinearArray}) = Broadcast.ArrayStyle{RectlinearArray}()

Base.similar(A::RectlinearArray{T,N}) where {T,N} = RectlinearArray(similar(A.data, A.dims), A.fixed_index)

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{RectlinearArray}}, ::Type{ElType}) where {ElType}
    A = find_ra(bc)
    RectlinearArray(similar(A.data, ElType, axes(bc)), A.fixed_index)
end

find_ra(bc::Base.Broadcast.Broadcasted) = find_ra(bc.args)
find_ra(args::Tuple) = find_ra(find_ra(args[1]), Base.tail(args))
find_ra(x) = x
find_ra(::Tuple{}) = nothing
find_ra(a::RectlinearArray, rest) = a
find_ra(::Any, rest) = find_ra(rest)

IndexStyle(::Type{<:RectlinearArray}) = IndexStyle(typeof(RectlinearArray).parameters[1])

# --- Begin Helper functions --- #

function _drop_index(indices::NTuple{N, Int}, k::Tuple) where {N}
    return Tuple([indices[i] for i in 1:N if i ∉ k])
end

function _skip_index(A::AbstractArray, fixed_index::Tuple, inds::Int...)
    skip = _drop_index(inds, fixed_index)
    return getindex(A, skip...)
end

function _drop_dims(A::AbstractArray, dims::Tuple)
    new_dims = ntuple(i -> i ∈ dims ? 1 : Colon(), Val(ndims(A)))
    return A[new_dims...]
end

# --- End Helper functions --- #

end
