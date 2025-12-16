# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

"""
    cell_center_gradient(mesh::SphericalGrid3D, A, I)

Compute the physical gradient ∇A at cell center `I` for scalar field `A`, returning an `SVector` of derivatives in the radial, polar, and azimuthal directions.
"""
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

"""
    edge_gradient(mesh::SphericalGrid3D, A, I)

Compute the physical gradient at the edge between cell `I` and its neighbor in each coordinate direction by averaging the corresponding edge derivatives.
"""
function edge_gradient(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}
) where {T}
  ∂A∂r = edge_derivative(mesh, A, I, 1)
  ∂A∂θ = edge_derivative(mesh, A, I, 2)
  ∂A∂ϕ = edge_derivative(mesh, A, I, 3)
  return @SVector [∂A∂r, ∂A∂θ, ∂A∂ϕ]
end
