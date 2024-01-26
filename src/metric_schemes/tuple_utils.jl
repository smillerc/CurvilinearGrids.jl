# Some simple functions to create the nested tuples that
# are used by the metric schemes

# Why not use (ξ, η, ζ) and (x, y, z) you may ask? Why
# make it more confusing? Becuase sometimes you may have
# a cylindrical 2D R-Z grid... and you want to access 
# R and Z which are x₁, x₂. This naming makes it flexible!
# And I'm sure you'll get over it! 😆
ξ_names(::NTuple{1,Int}) = (:∂ξ₁,)
ξ_names(::NTuple{2,Int}) = (:∂ξ₁, :∂ξ₂)
ξ_names(::NTuple{3,Int}) = (:∂ξ₁, :∂ξ₂, :∂ξ₃)

ξ̂_names(::NTuple{1,Int}) = (:∂ξ̂₁,)
ξ̂_names(::NTuple{2,Int}) = (:∂ξ̂₁, :∂ξ̂₂)
ξ̂_names(::NTuple{3,Int}) = (:∂ξ̂₁, :∂ξ̂₂, :∂ξ̂₃)

x_names(::NTuple{1,Int}) = (:∂x₁,)
x_names(::NTuple{2,Int}) = (:∂x₁, :∂x₂)
x_names(::NTuple{3,Int}) = (:∂x₁, :∂x₂, :∂x₃)

edge_names(::NTuple{1,Int}) = (:i₊½,)
edge_names(::NTuple{2,Int}) = (:i₊½, :j₊½)
edge_names(::NTuple{3,Int}) = (:i₊½, :j₊½, :k₊½)

function nested_tuple(dims, name_func, T=Float64)
  names = name_func(dims)
  return NamedTuple{(names)}(zeros(T, dims) for i in 1:length(names))
end