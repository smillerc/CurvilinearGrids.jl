module RectlinearArrays

using Adapt
using KernelAbstractions

export RectlinearArray

"""
    RectlinearArray{T,N,K} <: AbstractArray{T,N}

`N`-dimensional dense array with elements of type `T` and `K` fixed indices.

Indices that are fixed are not stored by the underlying data, providing performance and storage benefits over a traditional `N`-dimensional dense array, while still behaving like a standard `N`-dimensional array. Note: if a `RectlinearArray` is constructed in the REPL with underlying data being a GPU array, the REPL will not be able to display the array.
"""
struct RectlinearArray{T,N,D,A<:AbstractArray{T,D},K,M} <: AbstractArray{T,N}
    data::A
    dims::NTuple{N,Int}
    fixed_indices::NTuple{K,Int}
    valid_indices::NTuple{M,Int}
end

# --- Begin special constructor functions --- #

# Define a way to create an uninitalized RectlinearArray
"""
    RectlinearArray(data::AbstractArray{T}, fixed_indices::Tuple)

Construct a `RectlinearArray` from an array `A` with indices `fixed_indices` constant.

# Examples
```jldoctest
julia> A = [1 2; 1 2]
2×2 Matrix{Int64}:
 1  2
 1  2

julia> B = RectlinearArray(A, (1,))
2×2 RectlinearArray{Int64, 2, 1, Vector{Int64}, 1, 1}:
 1  2
 1  2
```
"""
function RectlinearArray(A::AbstractArray{T}, fixed_indices::Tuple) where {T}
    @assert ndims(A) > 0 "RectlinearArray does not support 0-dimensional arrays."
    new_A = _drop_dims(A, fixed_indices)
    valid_indices = filter(i -> i ∉ fixed_indices, Tuple(1:ndims(A)))
    return RectlinearArray{eltype(A), ndims(A), ndims(new_A), typeof(new_A), length(fixed_indices), length(valid_indices)}(new_A, size(A), fixed_indices, valid_indices)
end

"""
    RectlinearArray(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims...)
    RectlinearArray(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Tuple)

Construct an uninitalized `RectlinearArray` with type `T` and `fixed_indices`.

If the type is not specified, it defaults to `Float64`. If `backend` is not specified, it defaults to a typical CPU array.

# Examples
```jldoctest
julia> A = RectlinearArray(Float64, (1,), (2,2))
2×2 RectlinearArray{Float64, 2, 1, Vector{Float64}, 1, 1}:
 6.912e-310  6.912e-310
 6.912e-310  6.912e-310
```
"""
function RectlinearArray(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Dims) where {T}
    @assert length(dims) > 0 "RectlinearArray does not support 0-dimensional arrays."
    A = Array{T}(undef, dims)
    new_A = _drop_dims(A, fixed_indices)
    valid_indices = filter(i -> i ∉ fixed_indices, Tuple(1:ndims(A)))
    return RectlinearArray{eltype(A), ndims(A), ndims(new_A), typeof(new_A), length(fixed_indices), length(valid_indices)}(new_A, size(A), fixed_indices, valid_indices)
end
function RectlinearArray(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Int...) where {T}
    @assert length(dims) > 0 "RectlinearArray does not support 0-dimensional arrays."
    A = Array{T}(undef, dims)
    new_A = _drop_dims(A, fixed_indices)
    valid_indices = filter(i -> i ∉ fixed_indices, Tuple(1:ndims(A)))
    return RectlinearArray{eltype(A), ndims(A), ndims(new_A), typeof(new_A), length(fixed_indices), length(valid_indices)}(new_A, size(A), fixed_indices, valid_indices)
end
function RectlinearArray(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Dims) where {T}
    @assert length(dims) > 0 "RectlinearArray does not support 0-dimensional arrays."
    A = allocate(backend, T, dims)
    new_A = _drop_dims(A, fixed_indices)
    valid_indices = filter(i -> i ∉ fixed_indices, Tuple(1:ndims(A)))
    return RectlinearArray{eltype(A), ndims(A), ndims(new_A), typeof(new_A), length(fixed_indices), length(valid_indices)}(new_A, size(A), fixed_indices, valid_indices)
end
function RectlinearArray(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Int...) where {T}
    @assert length(dims) > 0 "RectlinearArray does not support 0-dimensional arrays."
    A = allocate(backend, T, dims)
    new_A = _drop_dims(A, fixed_indices)
    valid_indices = filter(i -> i ∉ fixed_indices, Tuple(1:ndims(A)))
    return RectlinearArray{eltype(A), ndims(A), ndims(new_A), typeof(new_A), length(fixed_indices), length(valid_indices)}(new_A, size(A), fixed_indices, valid_indices)
end
function RectlinearArray(fixed_indices::Tuple{Vararg{Int}}, dims::Dims)
    @assert length(dims) > 0 "RectlinearArray does not support 0-dimensional arrays."
    A = Array{Float64}(undef, dims)
    new_A = _drop_dims(A, fixed_indices)
    valid_indices = filter(i -> i ∉ fixed_indices, Tuple(1:ndims(A)))
    return RectlinearArray{eltype(A), ndims(A), ndims(new_A), typeof(new_A), length(fixed_indices), length(valid_indices)}(new_A, size(A), fixed_indices, valid_indices)
end
function RectlinearArray(fixed_indices::Tuple{Vararg{Int}}, dims::Int...)
    @assert length(dims) > 0 "RectlinearArray does not support 0-dimensional arrays."
    A = Array{Float64}(undef, dims)
    new_A = _drop_dims(A, fixed_indices)
    valid_indices = filter(i -> i ∉ fixed_indices, Tuple(1:ndims(A)))
    return RectlinearArray{eltype(A), ndims(A), ndims(new_A), typeof(new_A), length(fixed_indices), length(valid_indices)}(new_A, size(A), fixed_indices, valid_indices)
end

# Define a zeros function
function zeros(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Int...) where {T}
    return RectlinearArray(Base.zeros(T, dims...), fixed_indices)
end

function zeros(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Dims) where {T}
    return RectlinearArray(Base.zeros(T, dims), fixed_indices)
end
function zeros(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Int...) where {T}
    return RectlinearArray(KernelAbstractions.zeros(backend, T, dims...), fixed_indices)
end

function zeros(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Dims) where {T}
    return RectlinearArray(KernelAbstractions.zeros(backend, T, dims), fixed_indices)
end
function zeros(fixed_indices::Tuple{Vararg{Int}}, dims::Int...)
    return RectlinearArray(Base.zeros(Float64, dims...), fixed_indices)
end

function zeros(fixed_indices::Tuple{Vararg{Int}}, dims::Dims)
    return RectlinearArray(Base.zeros(Float64, dims), fixed_indices)
end

# --- End special constructor functions --- #

# --- Begin overload functions --- #

# Specify the size of the RectlinearArray
Base.size(A::RectlinearArray) = A.dims

# Define the getindex function
@inline function Base.getindex(A::RectlinearArray{T,N}, I::Vararg{Int, N}) where {T,N}
    _skip_index(A.data, A.valid_indices, I...)
end
# This getindex function is intended for indices that return multiple values (e.g. a range like 2:9). This allocates a new array (much like Julia's standard getindex in this case) and therefore will not work in GPU code
@inline function Base.getindex(A::RectlinearArray{T,N}, I::Vararg{Any, N}) where {T,N}
    return RectlinearArray(reshape([A[i...] for i in Iterators.product(I...)], length.(I)), A.fixed_indices)
end
# This getindex function is intended for indices that return multiple values (e.g. a CartesianIndices). This allocates a new array (much like Julia's standard getindex in this case) and therefore will not work in GPU code
@inline function Base.getindex(A::RectlinearArray{T,N}, I::CartesianIndices) where {T,N}
    return RectlinearArray(reshape([A[i] for i in I], size(I)), A.fixed_indices)
end

# Define the setindex! function
@inline function Base.setindex!(A::RectlinearArray{T,N}, v, I::Vararg{Int,N}) where {T,N}
    A.data[_drop_index(I, A.valid_indices)...] = v
end

# Specify the resultant array when broadcasting with .
Base.BroadcastStyle(::Type{<:RectlinearArray}) = Broadcast.ArrayStyle{RectlinearArray}()

# Define the similar function
Base.similar(A::RectlinearArray{T,N}) where {T,N} = RectlinearArray(similar(A.data, A.dims), A.fixed_indices)

# Specify how similar produces arrays for broadcasting
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

# Specify how to index a RectlinearArray
IndexStyle(::Type{<:RectlinearArray}) = IndexStyle(typeof(RectlinearArray).parameters[4])

# Overload eachindex for performance
Base.eachindex(A::RectlinearArray) = eachindex(A.data)

# Overload CartesianIndices for performance
# Base.CartesianIndices(A::RectlinearArray) = CartesianIndices(A.data)

# Overload materialize! for in-place broadcasting; Changed broadcasting. This now essentially operates over the "largest" eachindex for the given arguments
function Base.Broadcast.materialize!(dest::RectlinearArray, bc::Broadcast.Broadcasted)
    bc = Broadcast.instantiate(bc)
    ref = _largest_eachindex(dest, bc.args...)
    ref_index_space = eachindex(ref)

    for i in ref_index_space
        dest[i] = bc.f((arg isa AbstractArray ? arg[i] : arg for arg in bc.args)...)
    end
    return dest
end

# --- End overload functions --- #

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

function _largest_eachindex(A...)
    maxlen = -1
    best = nothing
    for a in A
        if isa(a, AbstractArray)
            len = length(eachindex(a))
            if len > maxlen
                maxlen = len
                best = a
            end
        end
    end
    return best
end

# --- End Helper functions --- #

Adapt.@adapt_structure RectlinearArray

end
