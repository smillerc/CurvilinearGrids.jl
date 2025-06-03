module SpecialArrays

using Adapt

export RectlinearArray

struct RectlinearArray{T,N,D,A<:AbstractArray{T,D},M} <: AbstractArray{T,N}
    data::A
    dims::NTuple{N,Int}
    fixed_index::NTuple{M,Int}
end

function RectlinearArray(data::AbstractArray{T}, fixed_index::Tuple) where {T}
    new_data = _drop_dims(data, fixed_index)
    new_data = isa(new_data, Number) ? [new_data] : new_data
    return RectlinearArray{eltype(data), ndims(data), ndims(new_data), typeof(new_data), length(fixed_index)}(new_data, size(data), fixed_index)
end

Base.size(A::RectlinearArray) = A.dims

@inline function Base.getindex(A::RectlinearArray{T,N}, I::Vararg{Int, N}) where {T,N}
    if ndims(A) == 1
        return A.data[1]
    elseif ndims(A) == 2
        if A.fixed_index == (1,)
            return A.data[I[2]]
        elseif A.fixed_index == (2,)
            return A.data[I[1]]
        else
            return A.data[1]
        end
    else
        if length(A.fixed_index) == 1
            if A.fixed_index == (1,)
                return A.data[I[2],I[3]]
            elseif A.fixed_index == (2,)
                return A.data[I[1],I[3]]
            else
                return A.data[I[1],I[2]]
            end
        elseif length(A.fixed_index) == 2
            if A.fixed_index == (1,2)
                return A.data[I[3]]
            elseif A.fixed_index == (2,3)
                return A.data[I[1]]
            else
                return A.data[I[2]]
            end
        else
            return A.data[1]
        end
    end
end

@inline function Base.setindex!(A::RectlinearArray{T,N}, v, I::Vararg{Int,N}) where {T,N}
    if ndims(A) == 1
        A.data[1] = v
    elseif ndims(A) == 2
        if A.fixed_index == (1,)
            A.data[I[2]] = v
        elseif A.fixed_index == (2,)
            A.data[I[1]] = v
        else
            A.data[1] = v
        end
    else
        if length(A.fixed_index) == 1
            if A.fixed_index == (1,)
                A.data[I[2],I[3]] = v
            elseif A.fixed_index == (2,)
                A.data[I[1],I[3]] = v
            else
                A.data[I[1],I[2]] = v
            end
        elseif length(A.fixed_index) == 2
            if A.fixed_index == (1,2)
                A.data[I[3]] = v
            elseif A.fixed_index == (2,3)
                A.data[I[1]] = v
            else
                A.data[I[2]] = v
            end
        else
            A.data[1] = v
        end
    end
end

Base.BroadcastStyle(::Type{<:RectlinearArray}) = Broadcast.ArrayStyle{RectlinearArray}()

Base.similar(A::RectlinearArray{T,N}) where {T,N} = RectlinearArray(similar(A.data, A.dims), A.fixed_index)

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{RectlinearArray}}, ::Type{ElType}) where {ElType}
    A = find_ra(bc)
    RectlinearArray(similar(A.data, A.dims), A.fixed_index)
end

find_ra(bc::Base.Broadcast.Broadcasted) = find_ra(bc.args)
find_ra(args::Tuple) = find_ra(find_ra(args[1]), Base.tail(args))
find_ra(x) = x
find_ra(::Tuple{}) = nothing
find_ra(a::RectlinearArray, rest) = a
find_ra(::Any, rest) = find_ra(rest)

IndexStyle(::Type{<:RectlinearArray}) = IndexStyle(typeof(RectlinearArray).parameters[1])

function Base.Broadcast.materialize!(dest::RectlinearArray, bc::Broadcast.Broadcasted)
    new_bc = Broadcast.instantiate(Broadcast.Broadcasted(bc.f, map(arg -> isa(arg, RectlinearArray) ? arg.data : arg, bc.args), axes(dest.data)))
    Broadcast.materialize!(dest.data, new_bc)
    return dest
end

# --- Begin Helper functions --- #
function _drop_dims(A::AbstractArray, dims::Tuple)
    new_dims = ntuple(i -> i ∈ dims ? 1 : Colon(), Val(ndims(A)))
    return A[new_dims...]
end
# --- End Helper functions --- #

Adapt.@adapt_structure RectlinearArray

end
