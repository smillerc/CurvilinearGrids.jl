# Functions used to work on CartesianIndex and CartesianIndices 
# to make looping and stencil creation easy for arbitary dimensions
# and axes

"""
  Apply a delta function to the cartesian index on a specified axis. For 
  example, `δ(3, CartesianIndex(1,2,3))` will give `CartesianIndex(0,0,1)`.
"""
δ(axis, ::CartesianIndex{N}) where {N} = CartesianIndex(ntuple(j -> j == axis ? 1 : 0, N))

"""
Get indices of of `±n` for an arbitary `axis`. For a 2d index, and
and offset on the `j` axis, this would be [i,j-n:j+n]. This makes
it simple to generate arbitrary stencils on arbitrary axes.

# Arguments
 - I::CartesianIndex{N}
 - axis::Int: which axis to provide the ±n to
 - n::Int: how much to offset on a given axis

 # Example
```julia
julia> I = CartesianIndex((4,5))
CartesianIndex(4, 5)

julia> plus_minus(I, 2, 2)
CartesianIndices((4:4, 3:7))
```

"""
function plus_minus(I::CartesianIndex{N}, axis::Int, n::Int) where {N}
  return (I - n * δ(axis, I)):(I + n * δ(axis, I))
end

"""
Go up by `n` on a given `axis` of a CartesianIndex
"""
function up(I::CartesianIndex{N}, axis::Int, n::Int) where {N}
  return I + n * δ(axis, I)
end

"""
Go down by `n` on a given `axis` of a CartesianIndex
"""
function down(I::CartesianIndex{N}, axis::Int, n::Int) where {N}
  return I - n * δ(axis, I)
end

"""
Expand the CartesianIndices ranges by `n` on `axis`
"""
function expand(
  domain::CartesianIndices{N,NTuple{N,UnitRange{T}}}, n::Int, axis::Int
) where {N,T}
  axis_mask = ntuple(j -> j == axis ? 1 : 0, N)
  return CartesianIndices(
    UnitRange.(first.(domain.indices) .- axis_mask, last.(domain.indices) .+ axis_mask)
  )
end

"""
Expand the CartesianIndices ranges by `n` on all axes
"""
function expand(domain::CartesianIndices{N,NTuple{N,UnitRange{T}}}, n::Int) where {N,T}
  return CartesianIndices(
    UnitRange.(first.(domain.indices) .- n, last.(domain.indices) .+ n)
  )
end

"""
Expand the CartesianIndices starting index by `-n` on all axes
"""
function expand_lower(
  domain::CartesianIndices{N,NTuple{N,UnitRange{T}}}, n::Int
) where {N,T}
  return CartesianIndices(UnitRange.(first.(domain.indices) .- n, last.(domain.indices)))
end

"""
Expand the CartesianIndices ending index by `+n` on all axes
"""
function expand_upper(
  domain::CartesianIndices{N,NTuple{N,UnitRange{T}}}, n::Int
) where {N,T}
  return CartesianIndices(UnitRange.(first.(domain.indices), last.(domain.indices) .+ n))
end