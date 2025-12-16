# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

"""
    cell_center_derivative(mesh::SphericalGrid3D, A, I, axis)

Return the physical derivative of the scalar field `A` at cell `I` along the `axis` direction. The derivative is computed from face-averaged values and effective distances derived from neighboring face areas and the local cell volume.
"""
function cell_center_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}
  I₊ = CartesianDomains.shift(I, axis, +1)
  I₋ = CartesianDomains.shift(I, axis, -1)

  @inbounds A₊ = A[I₊]
  @inbounds A₋ = A[I₋]

  Vᵢ = cell_volume(mesh, I)

  A_face_p = face_area_p(mesh, I, axis)   # A_{i+1/2}
  A_face_m = face_area_m(mesh, I, axis)   # A_{i-1/2}

  # effective distances using only current cell
  Δx_f = Vᵢ / A_face_p
  Δx_b = Vᵢ / A_face_m

  # second-order central derivative 
  ∂A = (A₊ - A₋) / (Δx_f + Δx_b)

  if axis == 3 # ϕ
    j = I.I[2]
    θ = mesh.centroid_coordinates.θ[j]
    ∂A = ∂A / sin(θ)
  end

  return ∂A
end

function cell_center_coordinate_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}
  I₊ = CartesianDomains.shift(I, axis, +1)
  I₋ = CartesianDomains.shift(I, axis, -1)

  @inbounds A₊ = A[I₊]
  @inbounds A₋ = A[I₋]

  Vᵢ = cell_volume(mesh, I)

  A_face_p = face_area_p(mesh, I, axis)   # A_{i+1/2}
  A_face_m = face_area_m(mesh, I, axis)   # A_{i-1/2}

  # Effective distance
  Δs = Vᵢ * (1 / A_face_p + 1 / A_face_m)

  # Second-order centered derivative
  return (A₊ - A₋) / Δs
end

# --------------------------------------------------------------------
# Edge operators
# --------------------------------------------------------------------

"""
    edge_derivative(mesh::SphericalGrid3D, A, I, axis)

Compute the derivative of `A` at the edge between cell `I` and its neighbor in `axis` by averaging the adjacent cell-centered derivatives. Useful for constructing edge-centered gradients and fluxes.
"""
function edge_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}

  #
  ∂Aᵢ = cell_center_derivative(mesh, A, I, axis)

  ᵢ₊₁ = CartesianDomains.shift(I, axis, +1)
  ∂Aᵢ₊₁ = cell_center_derivative(mesh, A, ᵢ₊₁, axis)

  ∂Aᵢ₊½ = (∂Aᵢ + ∂Aᵢ₊₁) / 2
  # 

  return ∂Aᵢ₊½
end
