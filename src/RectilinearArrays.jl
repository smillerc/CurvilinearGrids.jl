module RectilinearArrays

using Adapt
using KernelAbstractions
using GPUArrays

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

const AnyRectilinearArray{T,N} = Union{
  RectilinearArray{T,N,D,A,K,M} where {D,A<:AbstractArray{T,D},K,M},
  WrappedArray{
    T,N,RectilinearArray,RectilinearArray{T,N,D,A,K,M}
  } where {D,A<:AbstractArray{T,D},K,M},
  SubArray{
    T,N,RectilinearArray{T,N,D,A,K,M},<:Tuple{Vararg{Union{Integer,UnitRange{Int}},N}},false
  } where {D,A<:AbstractArray{T,D},K,M},
}

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
    T,ndims(A),ndims(new_A),typeof(new_A),K,ndims(A) - length(fixed_indices)
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
    T,
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(dims) - length(fixed_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end

function RectilinearArray(
  ::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Vararg{Int,N}
) where {T,N}
  @assert length(dims) > 0 "RectilinearArray does not support 0-dimensional arrays."
  A = Array{T}(undef, dims)
  new_A = _drop_dims(A, fixed_indices, CPU())
  valid_indices = filter(i -> i ∉ fixed_indices, ntuple(i -> i, N))
  return RectilinearArray{
    T,
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(dims) - length(fixed_indices),
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
    T,
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(dims) - length(fixed_indices),
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
    T,
    ndims(A),
    ndims(new_A),
    typeof(new_A),
    length(fixed_indices),
    length(dims) - length(fixed_indices),
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
    length(dims) - length(fixed_indices),
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
    length(dims) - length(fixed_indices),
  }(
    new_A, size(A), fixed_indices, valid_indices
  )
end

# Show function for the REPL
function Base.show(
  io::IO, ::MIME"text/plain", A::RectilinearArray{T,N,M,DA,Ax1,Ax2}
) where {T,N,M,DA<:AbstractGPUArray,Ax1,Ax2}
  summary_str = Base.summary(A)
  new_data = Array(A.data)
  A_show = RectilinearArray{T,N,ndims(new_data),typeof(new_data),Ax1,Ax2}(
    new_data, size(A), A.fixed_indices, A.valid_indices
  )
  io = IOContext(io, :compact => false)
  println(io, summary_str)
  Base.print_array(io, A_show)
end

function Base.show(
  io::IO, ::MIME"text/plain", A::SubArray{J,K,R,I,L}
) where {J,K,I,L,T,N,M,DA<:AbstractGPUArray,Ax1,Ax2,R<:RectilinearArray{T,N,M,DA,Ax1,Ax2}}
  A = parent(A)
  summary_str = Base.summary(A)
  new_data = Array(A.data)
  A_show = RectilinearArray{T,N,ndims(new_data),typeof(new_data),Ax1,Ax2}(
    new_data, size(A), A.fixed_indices, A.valid_indices
  )
  io = IOContext(io, :compact => false)
  println(io, summary_str)
  Base.print_array(io, A_show)
end

# Define a zeros function
@inline function zeros(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Int...) where {T}
  return RectilinearArray(Base.zeros(T, dims...), fixed_indices)
end

@inline function zeros(::Type{T}, fixed_indices::Tuple{Vararg{Int}}, dims::Dims) where {T}
  return RectilinearArray(Base.zeros(T, dims), fixed_indices)
end

@inline function zeros(
  ::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Int...
) where {T}
  return RectilinearArray(KernelAbstractions.zeros(backend, T, dims...), fixed_indices)
end

@inline function zeros(
  ::Type{T}, backend::Backend, fixed_indices::Tuple{Vararg{Int}}, dims::Dims
) where {T}
  return RectilinearArray(KernelAbstractions.zeros(backend, T, dims), fixed_indices)
end

@inline function zeros(fixed_indices::Tuple{Vararg{Int}}, dims::Int...)
  return RectilinearArray(Base.zeros(Float64, dims...), fixed_indices)
end

@inline function zeros(fixed_indices::Tuple{Vararg{Int}}, dims::Dims)
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

# getindex for GPU arrays. Scalar indexing is disallowed, but slicing is not
@inline function Base.getindex(
  A::RectilinearArray{T,N,D,DA,K,M}, inds::Vararg{Union{Colon,UnitRange{Int64}},N}
) where {T,N,D,DA<:AbstractGPUArray,K,M}
  formatted_inds = ntuple(
    i -> inds[i] isa Colon ? UnitRange(axes(A)[i]) : inds[i]::UnitRange{Int64}, Val(N)
  )
  gpu_view = view(A.data, _drop_index(formatted_inds, A.valid_indices)...)
  shape::NTuple{N,Base.OneTo{Int}} = Base.index_shape(formatted_inds...)
  result_size = ntuple(i -> shape[i].stop, Val(N))
  return RectilinearArray{T,N,ndims(gpu_view),typeof(parent(gpu_view)),K,M}(
    gpu_view, result_size, A.fixed_indices, A.valid_indices
  )
end

@inline function Base.getindex(
  A::RectilinearArray{T,N,D,DA,K,M}, inds::CartesianIndices
) where {T,N,D,DA<:AbstractGPUArray,K,M}
  is = ntuple(i -> i ∈ A.valid_indices ? UnitRange(axes(inds)[i]) : 1:1, ndims(A))
  vcor = view(inds, is...)
  vcors = ntuple(
    i -> vcor[1][A.valid_indices[i]]:vcor[end][A.valid_indices[i]], length(A.valid_indices)
  )
  gpu_view = view(A.data, vcors...)
  shape::NTuple{N,Base.OneTo{Int}} = Base.index_shape(inds)
  result_size = ntuple(i -> shape[i].stop, Val(N))
  return RectilinearArray{T,N,ndims(gpu_view),typeof(parent(gpu_view)),K,M}(
    gpu_view, result_size, A.fixed_indices, A.valid_indices
  )
end

# Define the setindex! function
@inline function Base.setindex!(A::RectilinearArray{T,N}, v, I::Vararg{Int,N}) where {T,N}
  A.data[_drop_index(I, A.valid_indices)...] = v
end

# KernelAbstractions backend function
function KernelAbstractions.get_backend(A::RectilinearArray)
  return KernelAbstractions.get_backend(A.data)
end

# This strides function will return the strides for a dense array with the size of A. 
Base.strides(A::RectilinearArray) = ntuple(i -> _prod_nminus1(size(A), i), length(size(A)))
function Base.eachindex(A::RectilinearArray)
  s = _drop_index(strides(A), A.valid_indices)
  a = _drop_index(axes(A), A.valid_indices)
  if isempty(a)
    return (1,)
  end
  return (sum(s[k] * (I[k] - 1) for k in eachindex(I)) + 1 for I in Iterators.product(a...))
end

function Base.CartesianIndices(A::RectilinearArray)
  CartesianIndices(_insert_ones_at(size(A.data), A.fixed_indices))
end

struct RectilinearArrayStyle{N} <: Broadcast.AbstractArrayStyle{N} end
RectilinearArrayStyle{N}(::Val{N}) where {N} = RectilinearArrayStyle{N}()

# Specify the resultant array when broadcasting with .
function Base.BroadcastStyle(::Type{<:RectilinearArray{T,N}}) where {T,N}
  RectilinearArrayStyle{N}()
end

function Base.BroadcastStyle(::Broadcast.ArrayStyle{<:RectilinearArray{T,N}}) where {T,N}
  RectilinearArrayStyle{N}()
end

function Base.BroadcastStyle(::RectilinearArrayStyle{N}, ::AbstractGPUArrayStyle) where {N}
  RectilinearArrayStyle{N}()
end

function Base.BroadcastStyle(
  ::Type{
    <:SubArray{
      T,
      N,
      RectilinearArray{T,N,D,DA,K,M},
      <:Tuple{Vararg{Union{Integer,UnitRange{Int}},N}},
      false,
    },
  },
) where {T,N,D,DA<:AbstractArray{T,D},K,M}
  RectilinearArrayStyle{N}()
end

function Base.BroadcastStyle(
  ::Broadcast.ArrayStyle{
    <:SubArray{
      T,
      N,
      RectilinearArray{T,N,D,DA,K,M},
      <:Tuple{Vararg{Union{Integer,UnitRange{Int}},N}},
      false,
    },
  },
) where {T,N,D,DA<:AbstractArray{T,D},K,M}
  RectilinearArrayStyle{N}()
end

# Define the similar function
function Base.similar(A::RectilinearArray{T,N}) where {T,N}
  RectilinearArray(similar(A.data, A.dims), A.fixed_indices)
end

# Specify how similar produces arrays for broadcasting
function Base.similar(
  bc::Broadcast.Broadcasted{RectilinearArrayStyle{N}}, ::Type{ElType}
) where {ElType,N}
  A = find_ra(bc)
  RectilinearArray(similar(A.data, A.dims), A.fixed_indices)
end

find_ra(bc::Base.Broadcast.Broadcasted) = find_ra(bc.args)
find_ra(args::Tuple) = find_ra(find_ra(args[1]), Base.tail(args))
find_ra(x) = x
find_ra(::Tuple{}) = nothing
find_ra(a::AnyRectilinearArray, rest) = parent(a)
find_ra(::Any, rest) = find_ra(rest)

arg_flatten(arg::AnyRectilinearArray, A) = parent(arg).data
arg_flatten(arg::Base.Broadcast.Broadcasted, A) = arg_flatten(Base.materialize(arg), A)
function arg_flatten(arg::AbstractArray, A)
  _drop_dims(arg, parent(A).fixed_indices, KernelAbstractions.get_backend(A))
end
arg_flatten(arg, A) = arg

function Base.broadcasted(
  f, bc::Broadcast.Broadcasted{RectilinearArrayStyle{N},T}
) where {T,N}
  A = find_ra(bc)
  Atype = typeof(A)
  data_args = ntuple(
    i -> begin
      if bc.args[i] isa RectilinearArray
        bc.args[i].data
      elseif bc.args[i] isa AbstractArray
        _drop_dims(bc.args[i], A.fixed_indices, KernelAbstractions.get_backend(A))
      else
        bc.args[i]
      end
    end,
    length(bc.args),
  )
  data_bc = Broadcast.materialize(Broadcast.broadcasted(f, data_args...))
  return RectilinearArray{Atype.parameters...}(
    data_bc, size(A), A.fixed_indices, A.valid_indices
  )
end

function Base.broadcast!(f, dest::RectilinearArray, args...)
  data_args = ntuple(i -> arg_flatten(args[i], dest), length(args))
  Base.broadcast!(f, dest.data, data_args...)
end

function copyto_assist(dest::SubArray)
  A = parent(dest)
  inds = dest.indices
  revised_inds = _drop_index(inds, A.valid_indices)
  return view(A.data, revised_inds...)
end
# function Base.copyto!(dest::SubArray{T,N,RectilinearArray{T,N,D,DA,K,M},<:Tuple{Vararg{Union{Integer,UnitRange{Int}},N}},false}, bc::Broadcast.Broadcasted{<:Broadcast.ArrayStyle{RectilinearArray}}) where {T,N,D,DA,K,M}
#     A = find_ra(bc)
# 
#     data_args = ntuple(i -> arg_flatten(bc.args[i], A), length(bc.args))
# 
#     raw_bc = Broadcast.broadcasted(bc.f, data_args...)
# 
#     if dest isa RectilinearArray
#         copyto!(dest.data, raw_bc)
#     else 
#         copyto!(copyto_assist(dest), raw_bc)
#     end
# 
#     return dest
# end

function Base.copyto!(
  dest::AnyRectilinearArray, bc::Broadcast.Broadcasted{RectilinearArrayStyle{N}}
) where {N}
  A = find_ra(bc)

  data_args = ntuple(i -> arg_flatten(bc.args[i], A), length(bc.args))

  raw_bc = Broadcast.broadcasted(bc.f, data_args...)

  if dest isa RectilinearArray
    copyto!(dest.data, raw_bc)
  else
    copyto!(copyto_assist(dest), raw_bc)
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

function _drop_index(indices::NTuple{N,UnitRange}, valid_indices::NTuple{M,Int}) where {N,M}
  return ntuple(i -> indices[valid_indices[i]], M)
end

function _drop_index(
  indices::NTuple{N,Base.OneTo{Int}}, valid_indices::NTuple{M,Int}
) where {N,M}
  return ntuple(i -> indices[valid_indices[i]], M)
end

function _skip_index(A::AbstractArray, valid_indices::Tuple, inds::Int...)
  skip = _drop_index(inds, valid_indices)
  return getindex(A, skip...)
end

@kernel function copy_kernel!(R, A, ::Val{dim_map}) where {dim_map}
  idx = @index(Global, Cartesian)
  full_idx = ntuple(i -> dim_map[i] > 0 ? idx[dim_map[i]] : 1, ndims(A))
  @inbounds R[idx] = A[full_idx...]
end

function _drop_dims(A::AbstractArray, dims::Tuple{Vararg{Int}}, backend::Backend)
  nd = ndims(A)
  all_dims = collect(1:nd)
  keep_dims = ntuple(j -> filter(i -> i ∉ dims, all_dims)[j], nd - length(dims))
  dim_map = ntuple(i -> i ∈ keep_dims ? findfirst(==(i), keep_dims) : 0, nd)

  new_shape = ntuple(i -> size(A, keep_dims[i]), nd - length(dims))
  R = similar(A, eltype(A), new_shape)

  kernel! = copy_kernel!(backend)
  kernel!(R, A, Val(dim_map); ndrange=size(R))

  return R
end

function _prod_nminus1(t::Tuple, n::Int)
  if n == 1
    return 1
  else
    acc = one(eltype(t))
    @inbounds for i in 1:(n - 1)
      acc *= t[i]
    end
  end
  return acc
end

function _insert_ones_at(
  I::Tuple{Vararg{Int,N}}, positions::Tuple{Vararg{Int,M}}
) where {N,M}
  return ntuple(j -> (j in positions ? 1 : I[j - count(p -> p < j, positions)]), Val(N + M))
end

# --- End Helper functions --- #

Adapt.@adapt_structure RectilinearArray

end