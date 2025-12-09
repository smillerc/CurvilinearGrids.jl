# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

function cell_center_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}
  I₊ = CartesianDomains.shift(I, axis, +1)
  I₋ = CartesianDomains.shift(I, axis, -1)

  A₊ = A[I₊]
  A₋ = A[I₋]

  # @show A₊, I₊
  # @show A₋, I₋

  Vᵢ = cell_volume(mesh, I)

  A_face_p = face_area_p(mesh, I, axis)   # A_{i+1/2}
  A_face_m = face_area_m(mesh, I, axis)   # A_{i-1/2}

  # @assert A_face_p > 0
  # @assert A_face_m > 0
  # @assert Vᵢ > 0

  # effective distances using only *current* volume
  Δx_f = Vᵢ / A_face_p
  Δx_b = Vᵢ / A_face_m

  # second-order central derivative using FV symmetric pair
  ∂A = (A₊ - A₋) / (Δx_f + Δx_b)

  if axis == 3 # ϕ
    j = I.I[2]
    θ = mesh.centroid_coordinates.θ[j]
    ∂A = ∂A / sin(θ)
  end

  return ∂A
end

# --------------------------------------------------------------------
# Edge operators
# --------------------------------------------------------------------

function edge_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int
) where {T}
  # I₊ = CartesianDomains.shift(I, axis, +1)
  # I₋ = CartesianDomains.shift(I, axis, -1)

  # @inbounds A₊ = A[I₊]
  # @inbounds A₋ = A[I₋]

  # Vᵢ = cell_volume(mesh, I)
  # A_edge = face_area_p(mesh, I, axis)   # A_{i+1/2}
  # Δs = Vᵢ / A_edge

  # ∂Aᵢ₊½ = (A₊ - A₋) / 2Δs

  # if axis == 3 # ϕ
  #   i, j, k = I.I
  #   # θ = mesh.node_coordinates.θ[j]
  #   θ = mesh.centroid_coordinates.θ[j]
  #   # r = mesh.centroid_coordinates.r[i]
  #   ∂Aᵢ₊½ = ∂Aᵢ₊½ / (sin(θ))
  # end
  if edge_on_lo_boundary(mesh, I, axis)
    ∂Aᵢ₊½ = one_sided_edge_derivative(mesh, A, I, axis, +2)
    # @show I, axis
    # error("bork")
  elseif edge_on_hi_boundary(mesh, I, axis)
    ∂Aᵢ₊½ = one_sided_edge_derivative(mesh, A, I, axis, -2)
  else
    ∂Aᵢ = cell_center_derivative(mesh, A, I, axis)

    ᵢ₊₁ = CartesianDomains.shift(I, axis, +1)
    ∂Aᵢ₊₁ = cell_center_derivative(mesh, A, ᵢ₊₁, axis)

    ∂Aᵢ₊½ = (∂Aᵢ + ∂Aᵢ₊₁) / 2
  end

  return ∂Aᵢ₊½
end

function one_sided_edge_derivative(
  mesh::SphericalGrid3D, A::AbstractArray{T,3}, I::CartesianIndex{3}, axis::Int, offset
) where {T}

  #  Identify first and second interior cells 
  I1 = CartesianDomains.shift(I, axis, sign(offset) * 1) # offset is ± 1 depending on the boundary
  I2 = CartesianDomains.shift(I, axis, sign(offset) * 2) # offset is ± 1 depending on the boundary

  #  cell-centered derivatives 
  # @show I1
  A1p = cell_center_derivative(mesh, A, I1, axis)
  # @show A1p
  # println()

  # @show I1
  A2p = cell_center_derivative(mesh, A, I2, axis)
  # @show A2p
  # println()

  coord = mesh.centroid_coordinates[axis]

  # --- Boundary face location ---
  node_coord = mesh.node_coordinates[axis]

  # Extract indices
  i, j, k = I.I
  i2, j2, k2 = I2.I

  #! format: off
  s1dim = axis == 1 ? i : axis == 2 ? j : k
  s2dim = axis == 1 ? i2 : axis == 2 ? j2 : k2
  #! format: on

  s_b = node_coord[s1dim]

  # --- Distances to interior cell centers ---
  @inbounds s1 = coord[s1dim]
  @inbounds s2 = coord[s2dim]

  h1 = abs(s1 - s_b)
  h2 = abs(s2 - s_b)

  # --- Quadratic reconstruction: 2nd-order boundary derivative ---
  return A1p + (A1p - A2p) * (h1 / (h1 + h2))
end
