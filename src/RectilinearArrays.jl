module RectilinearArrays

using Adapt
using KernelAbstractions

export RectilinearArray

"""
    RectilinearArray{T,N,K} <: AbstractArray{T,N}

`N`-dimensional dense array with elements of type `T` and `K` fixed indices.

Indices that are fixed are not stored by the underlying data, providing performance and storage benefits over a traditional `N`-dimensional dense array, while still behaving like a standard `N`-dimensional array. Note: if a `RectilinearArray` is constructed in the REPL with underlying data being a GPU array, the REPL will not be able to display the array.
"""
struct RectilinearArray{T,N,D,A<:AbstractArray{T,D},K,M} <: AbstractArray{T,N}
  data::A
  dims::NTuple{N,Int}
  fixed_indices::NTuple{K,Int}
  valid_indices::NTuple{M,Int}
end

# --- Begin special constructor functions --- #

# Define a way to create an uninitalized RectilinearArray
"""
    RectilinearArray(data::AbstractArray{T}, fixed_indices::Tuple)

Construct a `RectilinearArray` from an array `A` with indices `fixed_indices` constant.

# Examples
```jldoctest
julia> A = [1 2; 1 2]
2×2 Matrix{Int64}:
 1  2
 1  2

julia> B = RectilinearArray(A, (1,))
2×2 RectilinearArray{Int64, 2, 1, Vector{Int64}, 1, 1}:
 1  2
 1  2
```
"""
function RectilinearArray(A::AbstractArray{T}, fixed_indices::NTuple{K,Int}) where {T,K}
  @assert ndims(A) > 0 "RectilinearArray does not support 0-dimensional arrays."
  new_A = _drop_dims(A, fixed_indices, KernelAbstractions.get_backend(A))
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, ndims(A)))
  return RectilinearArray{
    eltype(A),
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    K,
    length(valid_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end

"""
    RectilinearArray(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims...)
    RectilinearArray(::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Tuple)

Construct an uninitalized `RectilinearArray` with type `T` and `fixed_indices`.

If the type is not specified, it defaults to `Float64`. If `backend` is not specified, it defaults to a typical CPU array.

# Examples
```jldoctest
julia> A = RectilinearArray(Float64, (1,), (2,2))
2×2 RectilinearArray{Float64, 2, 1, Vector{Float64}, 1, 1}:
 6.912e-310  6.912e-310
 6.912e-310  6.912e-310
```
"""
function RectilinearArray(
  ::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Dims{N}
) where {T,N}
  @assert length(dims) > 0 "RectilinearArray does not support 0-dimensional arrays."
  A = Array{T}(undef, dims)
  new_A = _drop_dims(A, fixed_indices, CPU())
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, N))
  return RectilinearArray{
    eltype(A),
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(valid_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end
function RectilinearArray(
    ::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Vararg{Int, N}
) where {T,N}
  @assert length(dims) > 0 "RectilinearArray does not support 0-dimensional arrays."
  A = Array{T}(undef, dims)
  new_A = _drop_dims(A, fixed_indices, CPU())
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, N))
  return RectilinearArray{
    eltype(A),
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(valid_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end
function RectilinearArray(
  ::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Dims{N}
) where {T,N}
  @assert length(dims) > 0 "RectilinearArray does not support 0-dimensional arrays."
  A = allocate(backend, T, dims)
  new_A = _drop_dims(A, fixed_indices, backend)
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, N))
  return RectilinearArray{
    eltype(A),
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(valid_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end
function RectilinearArray(
    ::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Vararg{Int,N}
) where {T,N}
  @assert length(dims) > 0 "RectilinearArray does not support 0-dimensional arrays."
  A = allocate(backend, T, dims)
  new_A = _drop_dims(A, fixed_indices, backend)
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, N))
  return RectilinearArray{
    eltype(A),
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(valid_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end
function RectilinearArray(fixed_indices::Tuple{Vararg{Int}}, dims::Dims{N}) where {N}
  @assert length(dims) > 0 "RectilinearArray does not support 0-dimensional arrays."
  A = Array{Float64}(undef, dims)
  new_A = _drop_dims(A, fixed_indices, CPU())
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, N))
  return RectilinearArray{
    eltype(A),
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(valid_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end
function RectilinearArray(fixed_indices::Tuple{Vararg{Int}}, dims::Vararg{Int,N}) where {N}
  @assert length(dims) > 0 "RectilinearArray does not support 0-dimensional arrays."
  A = Array{Float64}(undef, dims)
  new_A = _drop_dims(A, fixed_indices, CPU())
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, N))
  return RectilinearArray{
    eltype(A),
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(valid_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end

# Define a zeros function
function zeros(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Int...) where {T}
  return RectilinearArray(Base.zeros(T, dims...), fixed_indices)
end

function zeros(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Dims) where {T}
  return RectilinearArray(Base.zeros(T, dims), fixed_indices)
end
function zeros(
  ::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Int...
) where {T}
  return RectilinearArray(KernelAbstractions.zeros(backend, T, dims...), fixed_indices)
end

function zeros(
  ::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Dims
) where {T}
  return RectilinearArray(KernelAbstractions.zeros(backend, T, dims), fixed_indices)
end
function zeros(fixed_indices::Tuple{Vararg{Int}}, dims::Int...)
  return RectilinearArray(Base.zeros(Float64, dims...), fixed_indices)
end

function zeros(fixed_indices::Tuple{Vararg{Int}}, dims::Dims)
  return RectilinearArray(Base.zeros(Float64, dims), fixed_indices)
end

# --- End special constructor functions --- #

# --- Begin overload functions --- #

# Specify the size of the RectilinearArray
Base.size(A::RectilinearArray) = A.dims

# Define the getindex function
@inline function Base.getindex(A::RectilinearArray{T,N}, I::Vararg{Int,N}) where {T,N}
  _skip_index(A.data, A.valid_indices, I...)
end
# This getindex function is intended for indices that return multiple values (e.g. a range like 2:9). This allocates a new array (much like Julia's standard getindex in this case) and therefore will not work in GPU code
@inline function Base.getindex(
  A::RectilinearArray{T,N}, I::Vararg{UnitRange{Int},N}
) where {T,N}
  return RectilinearArray(
    reshape([A[i...] for i in Iterators.product(I...)], length.(I)), A.fixed_indices
  )
end
# This getindex function is intended for indices that return multiple values (e.g. a CartesianIndices). This allocates a new array (much like Julia's standard getindex in this case) and therefore will not work in GPU code
@inline function Base.getindex(A::RectilinearArray{T,N}, I::CartesianIndices) where {T,N}
  return RectilinearArray(reshape([A[i] for i in I], size(I)), A.fixed_indices)
end

# Define the setindex! function
@inline function Base.setindex!(A::RectilinearArray{T,N}, v, I::Vararg{Int,N}) where {T,N}
  A.data[_drop_index(I, A.valid_indices)...] = v
end

# Specify the resultant array when broadcasting with .
Base.BroadcastStyle(::Type{<:RectilinearArray}) = Broadcast.ArrayStyle{RectilinearArray}()

# Define the similar function
function Base.similar(A::RectilinearArray{T,N}) where {T,N}
  RectilinearArray(similar(A.data, A.dims), A.fixed_indices)
end

# Specify how similar produces arrays for broadcasting
function Base.similar(
  bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{RectilinearArray}}, ::Type{ElType}
) where {ElType}
  A = find_ra(bc)
  RectilinearArray(similar(A.data, A.dims), A.fixed_indices)
end

find_ra(bc::Base.Broadcast.Broadcasted) = find_ra(bc.args)
find_ra(args::Tuple) = find_ra(find_ra(args[1]), Base.tail(args))
find_ra(x) = x
find_ra(::Tuple{}) = nothing
find_ra(a::RectilinearArray, rest) = a
find_ra(::Any, rest) = find_ra(rest)

# Specify how to index a RectilinearArray
# Base.IndexStyle(::Type{<:RectilinearArray}) = IndexStyle(typeof(RectilinearArray).parameters[4])

# Overload eachindex for performance
Base.eachindex(A::RectilinearArray) = eachindex(A.data)

# Overload CartesianIndices for performance
# Base.CartesianIndices(A::RectilinearArray) = CartesianIndices(A.data)

# Overload materialize! for in-place broadcasting; Changed broadcasting. This now essentially operates over the "largest" eachindex for the given arguments
function Base.Broadcast.materialize!(dest::RectilinearArray, bc::Broadcast.Broadcasted)
  bc = Broadcast.instantiate(bc)
  ref = _largest_eachindex(dest, bc.args...)
  ref_index_space = eachindex(ref)

  for i in ref_index_space
    args_i = map(arg -> begin
      if arg isa AbstractArray
        arg[i]
      elseif arg isa Broadcast.Broadcasted
        Base.materialize(arg)[i]
      else
        arg
      end
    end, bc.args)
    dest[i] = bc.f(args_i...)
  end
  return dest
end

# Overload the all function to enforce boolean values
function Base.all(A::RectilinearArray)
    return all(Bool.(A.data))
end

# --- End overload functions --- #

# --- Begin Helper functions --- #

function _drop_index(indices::NTuple{N,Int}, valid_indices::NTuple{M,Int}) where {N,M}
  return ntuple(i -> indices[valid_indices[i]], M)
end

function _skip_index(A::AbstractArray, valid_indices::Tuple, inds::Int...)
  skip = _drop_index(inds, valid_indices)
  return getindex(A, skip...)
end

# function _drop_dims(A::AbstractArray, dims::Tuple)
#   new_dims = ntuple(i -> i ∈ dims ? 1 : Colon(), Val(ndims(A)))
#   return isa(A[new_dims...], Number) ? [A[new_dims...]] : A[new_dims...]
# end
@kernel function copy_kernel!(R, A, keep_dims)
    idx = @index(Global, Cartesian)
    full_idx = ntuple(i -> (i ∈ keep_dims ? idx[findfirst(==(i), keep_dims)] : 1), ndims(A))
    @inbounds R[idx] = A[full_idx...]
end
function _drop_dims(A::AbstractArray, dims::Tuple{Vararg{Int}}, backend::Backend)
    nd = ndims(A)
    all_dims = collect(1:nd)
    keep_dims = ntuple(j -> filter(i -> i ∉ dims, all_dims)[j], nd - length(dims))

    new_shape = ntuple(i -> size(A, keep_dims[i]), nd - length(dims))
    R = similar(A, eltype(A), new_shape)

    kernel! = copy_kernel!(backend)
    kernel!(R, A, keep_dims; ndrange=size(R))

    return R
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

Adapt.@adapt_structure RectilinearArray

end
