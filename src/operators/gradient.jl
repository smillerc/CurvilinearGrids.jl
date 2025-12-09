# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

# ∇A at the cell center (physical gradient)
function cell_center_gradient(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∂A∂r = cell_center_derivative(mesh, A, I, 1)
  ∂A∂θ = cell_center_derivative(mesh, A, I, 2)
  ∂A∂ϕ = cell_center_derivative(mesh, A, I, 3)
  return @SVector [∂A∂r, ∂A∂θ, ∂A∂ϕ]
end

# --------------------------------------------------------------------
# Edge operators
# --------------------------------------------------------------------

# ∇A at the edge center (full physical gradient)
function edge_gradient(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∂A∂r = edge_derivative(mesh, A, I, 1)
  ∂A∂θ = edge_derivative(mesh, A, I, 2)
  ∂A∂ϕ = edge_derivative(mesh, A, I, 3)
  return @SVector [∂A∂r, ∂A∂θ, ∂A∂ϕ]
end
