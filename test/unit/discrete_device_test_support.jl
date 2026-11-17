using CurvilinearGrids, KernelAbstractions, Adapt

# Device-shaped storage rejects scalar access, so host interpolation cannot
# accidentally index device buffers during construction or refresh.
struct DiscreteTestBackend <: KernelAbstractions.Backend end
struct DiscreteTestArray{T,N} <: AbstractArray{T,N}
  data::Array{T,N}
end
Base.size(a::DiscreteTestArray) = size(a.data)
Base.getindex(::DiscreteTestArray, I...) = error("device scalar read")
Base.setindex!(::DiscreteTestArray, value, I...) = error("device scalar write")
Base.Array(a::DiscreteTestArray) = copy(a.data)
Base.copyto!(a::DiscreteTestArray, b::Array) = (copyto!(a.data, b); a)
KernelAbstractions.get_backend(::DiscreteTestArray) = DiscreteTestBackend()
KernelAbstractions.allocate(::DiscreteTestBackend, ::Type{T}, dims::Tuple) where T = DiscreteTestArray(Array{T}(undef,dims))
KernelAbstractions.synchronize(::DiscreteTestBackend) = nothing
Adapt.adapt_storage(::Type{Array}, a::DiscreteTestArray) = Array(a)
struct DiscreteTestAdaptor end
Adapt.adapt_storage(::DiscreteTestAdaptor, a::Array) = DiscreteTestArray(copy(a))

host_data(a::AbstractArray) = Array(a)
host_data(a::Union{Tuple,NamedTuple}) = map(host_data,a)
