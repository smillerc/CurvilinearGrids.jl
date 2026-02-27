
# --------------------------------------------------------------------
# Cell-center operators
# --------------------------------------------------------------------

# ∇×A at the cell center (physical components)
function cell_center_curl(mesh::SphericalGrid3D, (A_r, A_θ, A_ϕ), I::CartesianIndex{3})
  @inbounds begin
    dAr_dθ = cell_center_derivative(mesh, A_r, I, 2)
    dAr_dϕ = cell_center_derivative(mesh, A_r, I, 3)
    dAθ_dϕ = cell_center_derivative(mesh, A_θ, I, 3)

    # ============================================================
    #    Composite derivatives:
    #    ∂θ(Aϕ sinθ),  ∂r(r Aϕ),  ∂r(r Aθ)
    # ============================================================

    Vᵢ = cell_volume(mesh, I)

    # ---------- ∂θ (Aϕ sinθ) ----------
    Iθ₊ = CartesianDomains.shift(I, 2, +1)
    Iθ₋ = CartesianDomains.shift(I, 2, -1)

    θ₊ = mesh.centroid_coordinates[2][Iθ₊.I[2]]
    θ₋ = mesh.centroid_coordinates[2][Iθ₋.I[2]]

    Aϕ_sinθ₊ = A_ϕ[Iθ₊] * sin(θ₊)
    Aϕ_sinθ₋ = A_ϕ[Iθ₋] * sin(θ₋)

    A_face_θ₊ = face_area_p(mesh, I, 2)
    A_face_θ₋ = face_area_m(mesh, I, 2)

    Δs_θ₊ = Vᵢ / A_face_θ₊
    Δs_θ₋ = Vᵢ / A_face_θ₋
    Δs_θ = Δs_θ₊ + Δs_θ₋

    dAϕsinθ_dθ = (Aϕ_sinθ₊ - Aϕ_sinθ₋) / Δs_θ

    # ---------- ∂r (r Aϕ) ----------
    Ir₊ = CartesianDomains.shift(I, 1, +1)
    Ir₋ = CartesianDomains.shift(I, 1, -1)

    r₊ = mesh.centroid_coordinates[1][Ir₊.I[1]]
    r₋ = mesh.centroid_coordinates[1][Ir₋.I[1]]

    A_face_r₊ = face_area_p(mesh, I, 1)
    A_face_r₋ = face_area_m(mesh, I, 1)

    Δs_r₊ = Vᵢ / A_face_r₊
    Δs_r₋ = Vᵢ / A_face_r₋
    Δs_r = Δs_r₊ + Δs_r₋

    rA_ϕ₊ = r₊ * A_ϕ[Ir₊]
    rA_ϕ₋ = r₋ * A_ϕ[Ir₋]

    drAϕ_dr = (rA_ϕ₊ - rA_ϕ₋) / Δs_r

    # ---------- ∂r (r Aθ) ----------
    rA_θ₊ = r₊ * A_θ[Ir₊]
    rA_θ₋ = r₋ * A_θ[Ir₋]

    drAθ_dr = (rA_θ₊ - rA_θ₋) / Δs_r

    # ============================================================
    # 3. Assemble spherical curl using physical derivatives
    # ============================================================
    (i, j, k) = I.I
    r = mesh.centroid_coordinates[1][i]
    θ = mesh.centroid_coordinates[2][j]
    sinθ = sin(θ)

    curl_r = (1 / (r * sinθ)) * (dAϕsinθ_dθ - dAθ_dϕ)
    curl_θ = (1 / r) * ((1 / sinθ) * dAr_dϕ - drAϕ_dr)
    curl_ϕ = (1 / r) * (drAθ_dr - dAr_dθ)

    return @SVector[curl_r, curl_θ, curl_ϕ]
  end
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
