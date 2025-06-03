module SpecialArrays

export RectlinearArray

mutable struct RectlinearArray{T,N,D,A<:Union{AbstractArray{T,D},T},M} <: AbstractArray{T,N}
    data::A
    dims::NTuple{N,Int}
    fixed_index::NTuple{M,Int}
end

function RectlinearArray(data::AbstractArray{T}, fixed_index::Tuple) where {T}
    new_data = _drop_dims(data, fixed_index)
    return RectlinearArray{eltype(data), ndims(data), ndims(new_data), typeof(new_data), length(fixed_index)}(new_data, size(data), fixed_index)
end

Base.size(A::RectlinearArray) = A.dims

@inline function Base.getindex(A::RectlinearArray{T,N}, I::Vararg{Int, N}) where {T,N}
    if !isa(A.data, Number)
        _skip_index(A.data, A.fixed_index, I...)
    else
        A.data
    end
end

@inline function Base.setindex!(A::RectlinearArray{T,N}, v, I::Vararg{Int,N}) where {T,N}
    if !isa(A.data, Number)
        A.data[_drop_index(I, A.fixed_index)...] = v
    else
        A.data = v
    end
end

Base.BroadcastStyle(::Type{<:RectlinearArray}) = Broadcast.ArrayStyle{RectlinearArray}()

Base.similar(A::RectlinearArray{T,N}) where {T,N} = RectlinearArray(similar((isa(A.data, Number) ? [A.data] : A.data), A.dims), A.fixed_index)

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{RectlinearArray}}, ::Type{ElType}) where {ElType}
    A = find_ra(bc)
    RectlinearArray(similar((isa(A.data, Number) ? [A.data] : A.data), A.dims), A.fixed_index)
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
