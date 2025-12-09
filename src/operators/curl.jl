
# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

# ∇×A at the cell center (physical components)
function cell_center_curl(
  mesh::SphericalGrid3D, (A_r, A_θ, A_ϕ)::NTuple{3,AbstractArray{T,3}}, I::CartesianIndex{3}
) where {T}
  V = cell_volume(mesh, I)

  I_i₊ = CartesianDomains.shift(I, 1, +1)
  I_i₋ = CartesianDomains.shift(I, 1, -1)
  I_j₊ = CartesianDomains.shift(I, 2, +1)
  I_j₋ = CartesianDomains.shift(I, 2, -1)
  I_k₊ = CartesianDomains.shift(I, 3, +1)
  I_k₋ = CartesianDomains.shift(I, 3, -1)

  # ω_r: circulation in (θ,ϕ)
  A_ϕ_face_θp = face_val(A_ϕ, I, I_j₊)
  A_ϕ_face_θm = face_val(A_ϕ, I_j₋, I)

  A_θ_face_ϕp = face_val(A_θ, I, I_k₊)
  A_θ_face_ϕm = face_val(A_θ, I_k₋, I)

  A_θ_area_p = face_area_p(mesh, I, 2)
  A_θ_area_m = face_area_m(mesh, I, 2)
  A_ϕ_area_p = face_area_p(mesh, I, 3)
  A_ϕ_area_m = face_area_m(mesh, I, 3)

  circ_r =
    (A_θ_area_p * A_ϕ_face_θp - A_θ_area_m * A_ϕ_face_θm) -
    (A_ϕ_area_p * A_θ_face_ϕp - A_ϕ_area_m * A_θ_face_ϕm)

  ω_r = circ_r / V

  # ω_θ: circulation in (ϕ,r)
  A_r_face_ϕp = face_val(A_r, I, I_k₊)
  A_r_face_ϕm = face_val(A_r, I_k₋, I)

  A_ϕ_face_rp = face_val(A_ϕ, I, I_i₊)
  A_ϕ_face_rm = face_val(A_ϕ, I_i₋, I)

  A_ϕ_area_p2 = face_area_p(mesh, I, 3)
  A_ϕ_area_m2 = face_area_m(mesh, I, 3)
  A_r_area_p = face_area_p(mesh, I, 1)
  A_r_area_m = face_area_m(mesh, I, 1)

  circ_θ =
    (A_ϕ_area_p2 * A_r_face_ϕp - A_ϕ_area_m2 * A_r_face_ϕm) -
    (A_r_area_p * A_ϕ_face_rp - A_r_area_m * A_ϕ_face_rm)

  ω_θ = circ_θ / V

  # ω_ϕ: circulation in (r,θ)
  A_θ_face_rp2 = face_val(A_θ, I, I_i₊)
  A_θ_face_rm2 = face_val(A_θ, I_i₋, I)

  A_r_face_θp2 = face_val(A_r, I, I_j₊)
  A_r_face_θm2 = face_val(A_r, I_j₋, I)

  A_r_area_p2 = face_area_p(mesh, I, 1)
  A_r_area_m2 = face_area_m(mesh, I, 1)
  A_θ_area_p2 = face_area_p(mesh, I, 2)
  A_θ_area_m2 = face_area_m(mesh, I, 2)

  circ_ϕ =
    (A_r_area_p2 * A_θ_face_rp2 - A_r_area_m2 * A_θ_face_rm2) -
    (A_θ_area_p2 * A_r_face_θp2 - A_θ_area_m2 * A_r_face_θm2)

  ω_ϕ = circ_ϕ / V

  return @SVector [ω_r, ω_θ, ω_ϕ]
end

# --------------------------------------------------------------------
# Edge operators
# --------------------------------------------------------------------

# ∇×A at the edge center = average of adjacent cell-center curls
function edge_curl(
  mesh::SphericalGrid3D,
  (A_r_edge, A_θ_edge, A_ϕ_edge)::NTuple{3,AbstractArray{T,3}},
  I::CartesianIndex{3},
  edge_axis::Int,

  
) where {T}
  I₊ = CartesianDomains.shift(I, edge_axis, +1)

  curl_L = cell_center_curl(mesh, (A_r_edge, A_θ_edge, A_ϕ_edge), I)
  curl_R = cell_center_curl(mesh, (A_r_edge, A_θ_edge, A_ϕ_edge), I₊)

  return (1 / 2) * (curl_L + curl_R)
end