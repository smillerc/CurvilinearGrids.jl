module SpecialArrays

export RectlinearArray

struct RectlinearArray{T,N,D,A<:AbstractArray{T,D},K,M} <: AbstractArray{T,N}
    data::A
    dims::NTuple{N,Int}
    fixed_indices::NTuple{K,Int}
    valid_indices::NTuple{M,Int}
end

function RectlinearArray(data::AbstractArray{T}, fixed_index::Tuple) where {T}
    new_data = _drop_dims(data, fixed_index)
    valid_indices = filter(i -> i ∉ fixed_index, Tuple(1:ndims(data)))
    return RectlinearArray{eltype(data), ndims(data), ndims(new_data), typeof(new_data), length(fixed_index), length(valid_indices)}(new_data, size(data), fixed_index, valid_indices)
end

Base.size(A::RectlinearArray) = A.dims

@inline function Base.getindex(A::RectlinearArray{T,N}, I::Vararg{Int, N}) where {T,N}
    _skip_index(A.data, A.valid_indices, I...)
end

@inline function Base.setindex!(A::RectlinearArray{T,N}, v, I::Vararg{Int,N}) where {T,N}
    A.data[_drop_index(I, A.valid_indices)...] = v
end

Base.BroadcastStyle(::Type{<:RectlinearArray}) = Broadcast.ArrayStyle{RectlinearArray}()

Base.similar(A::RectlinearArray{T,N}) where {T,N} = RectlinearArray(similar(A.data, A.dims), A.fixed_indices)

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{RectlinearArray}}, ::Type{ElType}) where {ElType}
    A = find_ra(bc)
    RectlinearArray(similar(A.data, A.dims), A.fixed_indices)
end

find_ra(bc::Base.Broadcast.Broadcasted) = find_ra(bc.args)
find_ra(args::Tuple) = find_ra(find_ra(args[1]), Base.tail(args))
find_ra(x) = x
find_ra(::Tuple{}) = nothing
find_ra(a::RectlinearArray, rest) = a
find_ra(::Any, rest) = find_ra(rest)

IndexStyle(::Type{<:RectlinearArray}) = IndexStyle(typeof(RectlinearArray).parameters[1])

function Base.Broadcast.materialize!(dest::RectlinearArray, bc::Broadcast.Broadcasted)
    if isa(dest.data, Number)
        dest.data = bc.f(map(arg -> isa(arg, RectlinearArray) ? arg.data : arg, bc.args)...)
    else
        new_bc = Broadcast.instantiate(Broadcast.Broadcasted(bc.f, map(arg -> isa(arg, RectlinearArray) ? arg.data : arg, bc.args), axes(dest.data)))
        Broadcast.materialize!(dest.data, new_bc)
    end 
    return dest
end

# --- Begin Helper functions --- #

function _drop_index(indices::NTuple{N, Int}, valid_indices::NTuple{M,Int}) where {N,M}
    return ntuple(i -> indices[valid_indices[i]], M)
end

function _skip_index(A::AbstractArray, valid_indices::Tuple, inds::Int...)
    skip = _drop_index(inds, valid_indices)
    return getindex(A, skip...)
end

function _drop_dims(A::AbstractArray, dims::Tuple)
    new_dims = ntuple(i -> i ∈ dims ? 1 : Colon(), Val(ndims(A)))
    return isa(A[new_dims...], Number) ? [A[new_dims...]] : A[new_dims...]
end

# --- End Helper functions --- #

end
