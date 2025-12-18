
# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

"""
    cell_center_divergence(mesh::SphericalGrid3D, (A_r, A_θ, A_ϕ), I)

Evaluate ∇⋅A at the cell center `I` for a vector field with physical components `(A_r, A_θ, A_ϕ)`. The divergence is formed from face-averaged fluxes scaled by the appropriate face areas and normalized by the cell volume.
"""
function cell_center_divergence(
  mesh::SphericalGrid3D, (A_r, A_θ, A_ϕ)::NTuple{3,AbstractArray{T,3}}, I::CartesianIndex{3}
) where {T}
  V = cell_volume(mesh, I)

  # convenience neighbor indices
  I_i₊ = CartesianDomains.shift(I, 1, +1)
  I_i₋ = CartesianDomains.shift(I, 1, -1)
  I_j₊ = CartesianDomains.shift(I, 2, +1)
  I_j₋ = CartesianDomains.shift(I, 2, -1)
  I_k₊ = CartesianDomains.shift(I, 3, +1)
  I_k₋ = CartesianDomains.shift(I, 3, -1)

  # radial faces
  @inbounds A_rᵢ = A_r[I]
  @inbounds A_rᵢ₊₁ = A_r[I_i₊]
  @inbounds A_rᵢ₋₁ = A_r[I_i₋]

  A_r_face_p = (1 / 2) * (A_rᵢ + A_rᵢ₊₁)
  A_r_face_m = (1 / 2) * (A_rᵢ₋₁ + A_rᵢ)

  A_r_area_p = face_area_p(mesh, I, 1)
  A_r_area_m = face_area_m(mesh, I, 1)

  flux_r = A_r_face_p * A_r_area_p - A_r_face_m * A_r_area_m

  # theta faces
  @inbounds A_θⱼ = A_θ[I]
  @inbounds A_θⱼ₊₁ = A_θ[I_j₊]
  @inbounds A_θⱼ₋₁ = A_θ[I_j₋]

  A_θ_face_p = (1 / 2) * (A_θⱼ + A_θⱼ₊₁)
  A_θ_face_m = (1 / 2) * (A_θⱼ₋₁ + A_θⱼ)

  A_θ_area_p = face_area_p(mesh, I, 2)
  A_θ_area_m = face_area_m(mesh, I, 2)

  flux_θ = A_θ_face_p * A_θ_area_p - A_θ_face_m * A_θ_area_m

  # phi faces
  @inbounds A_ϕₖ = A_ϕ[I]
  @inbounds A_ϕₖ₊₁ = A_ϕ[I_k₊]
  @inbounds A_ϕₖ₋₁ = A_ϕ[I_k₋]

  A_ϕ_face_p = (1 / 2) * (A_ϕₖ + A_ϕₖ₊₁)
  A_ϕ_face_m = (1 / 2) * (A_ϕₖ₋₁ + A_ϕₖ)

  A_ϕ_area_p = face_area_p(mesh, I, 3)
  A_ϕ_area_m = face_area_m(mesh, I, 3)

  flux_ϕ = A_ϕ_face_p * A_ϕ_area_p - A_ϕ_face_m * A_ϕ_area_m

  div_A = (flux_r + flux_θ + flux_ϕ) / V
  return div_A
end

# --------------------------------------------------------------------
# Edge operators
# --------------------------------------------------------------------

"""
    edge_divergence(mesh::SphericalGrid3D, (A_r, A_θ, A_ϕ), I, edge_axis)

Return the edge-centered divergence by averaging the cell-centered divergences on either side of the edge aligned with `edge_axis`.
"""
function edge_divergence(
  mesh::SphericalGrid3D,
  (A_r, A_θ, A_ϕ)::NTuple{3,AbstractArray{T,3}},
  I::CartesianIndex{3},
  edge_axis::Int,
) where {T}
  I₊ = CartesianDomains.shift(I, edge_axis, +1)

  div_L = cell_center_divergence(mesh, (A_r, A_θ, A_ϕ), I)
  div_R = cell_center_divergence(mesh, (A_r, A_θ, A_ϕ), I₊)

  return (1 / 2) * (div_L + div_R)
end
