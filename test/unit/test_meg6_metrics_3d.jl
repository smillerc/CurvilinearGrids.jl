using Test

using CurvilinearGrids
using StaticArrays
using Test
using BenchmarkTools
using WriteVTK

function save_vtk(mesh, I)
  fn = "wavy"
  @info "Writing to $fn.vti"

  xyz_n = CurvilinearGrids.coords(mesh)
  domain = mesh.iterators.cell.domain

  @views begin
    #   J = [m.J for m in mesh.cell_center_metrics[domain]]
    ξx = [m.x for m in mesh.cell_center_metrics.ξ[domain]]
    ξ̂xᵢ₊½ = [m.x for m in mesh.edge_metrics.i₊½.ξ̂[domain]]
    #   ηy = [m.ηy for m in mesh.cell_center_metrics[domain]]
  end

  vtk_grid(fn, xyz_n) do vtk

    # w1 = [x[1] for idx in eachindex(f.char_state.W[ilo:ihi, jlo:jhi]

    vtk["J"] = mesh.cell_center_metrics.J[domain]
    vtk["ξx"] = ξx
    vtk["ξ̂xᵢ₊½"] = ξ̂xᵢ₊½
    vtk["I"] = I[domain]
    # vtk["ηy"] = ηy
  end
end

# @testset "3D wavy mesh" 
# begin
function wavy_grid(ni, nj, nk)
  Lx = Ly = Lz = 4.0

  xmin = -Lx / 2
  ymin = -Ly / 2
  zmin = -Lz / 2

  Δx0 = Lx / ni
  Δy0 = Ly / nj
  Δz0 = Lz / nk

  x(i, j, k) = xmin + Δx0 * ((i - 1) + sinpi((j - 1) * Δy0) * sinpi((k - 1) * Δz0))
  y(i, j, k) = ymin + Δy0 * ((j - 1) + sinpi((k - 1) * Δz0) * sinpi((i - 1) * Δx0))
  z(i, j, k) = zmin + Δz0 * ((k - 1) + sinpi((i - 1) * Δx0) * sinpi((j - 1) * Δy0))

  return (x, y, z)
end

ni = nj = nk = 20
nhalo = 4
x, y, z = wavy_grid(ni, nj, nk)
mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)

domain = mesh.iterators.cell.domain

I₁_passes = true
I₂_passes = true
I₃_passes = true
I = similar(mesh.cell_center_metrics.J)
for idx in domain
  i, j, k = idx.I .+ 0.5
  # m_i₊½ = metrics(mesh, (i + 0.5, j, k), 0)
  # m_i₋½ = metrics(mesh, (i - 0.5, j, k), 0)
  # m_j₊½ = metrics(mesh, (i, j + 0.5, k), 0)
  # m_j₋½ = metrics(mesh, (i, j - 0.5, k), 0)
  # m_k₊½ = metrics(mesh, (i, j, k + 0.5), 0)
  # m_k₋½ = metrics(mesh, (i, j, k - 0.5), 0)

  # I₁ = (m_i₊½.ξ.x - m_i₋½.ξ.x) + (m_j₊½.η.x - m_j₋½.η.x) + (m_k₊½.ζ.x - m_k₋½.ζ.x)
  # I₂ = (m_i₊½.ξ.y - m_i₋½.ξ.y) + (m_j₊½.η.y - m_j₋½.η.y) + (m_k₊½.ζ.y - m_k₋½.ζ.y)
  # I₃ = (m_i₊½.ξ.z - m_i₋½.ξ.z) + (m_j₊½.η.z - m_j₋½.η.z) + (m_k₊½.ζ.z - m_k₋½.ζ.z)
  m_i₊½ = conservative_metrics(mesh, (i + 0.5, j, k))
  m_i₋½ = conservative_metrics(mesh, (i - 0.5, j, k))
  m_j₊½ = conservative_metrics(mesh, (i, j + 0.5, k))
  m_j₋½ = conservative_metrics(mesh, (i, j - 0.5, k))
  m_k₊½ = conservative_metrics(mesh, (i, j, k + 0.5))
  m_k₋½ = conservative_metrics(mesh, (i, j, k - 0.5))

  I₁ = (m_i₊½.ξ̂.x - m_i₋½.ξ̂.x) + (m_j₊½.η̂.x - m_j₋½.η̂.x) + (m_k₊½.ζ̂.x - m_k₋½.ζ̂.x)
  I₂ = (m_i₊½.ξ̂.y - m_i₋½.ξ̂.y) + (m_j₊½.η̂.y - m_j₋½.η̂.y) + (m_k₊½.ζ̂.y - m_k₋½.ζ̂.y)
  I₃ = (m_i₊½.ξ̂.z - m_i₋½.ξ̂.z) + (m_j₊½.η̂.z - m_j₋½.η̂.z) + (m_k₊½.ζ̂.z - m_k₋½.ζ̂.z)

  I[idx] = I₁
  I₁ = I₁ * (abs(I₁) > eps())
  I₂ = I₂ * (abs(I₂) > eps())
  I₃ = I₃ * (abs(I₃) > eps())

  # @show m_i₊½
  if I₁ > 0
    @show I₁, I₂, I₃
    @show i, j, k
    # @show m_i₊½.ξ̂.x, m_i₋½.ξ̂.x
    # @show m_j₊½.η̂.x, m_j₋½.η̂.x
    # @show m_k₊½.ζ̂.x, m_k₋½.ζ̂.x
    error("done")
  end

  I₁_passes = abs(I₁) < eps()
  I₂_passes = abs(I₂) < eps()
  I₃_passes = abs(I₃) < eps()
  # if !(I₁_passes && I₂_passes && I₃_passes)
  #   @warn "fail!"
  #   break
  # end
end
@test I₁_passes
@test I₂_passes
@test I₃_passes
# ∇y((i, j, k)) = ForwardDiff.gradient(x -> mesh.y_func(i, j, k), @SVector [i, j, k])
# ∇y((4, 5, 60))

# _x((i, j, k)) = mesh.x_func(i, j, k)
# _y((i, j, k)) = mesh.y_func(i, j, k)
# _z((i, j, k)) = mesh.z_func(i, j, k)

# # gradx_y((i, j, k)) = ForwardDiff.gradient(_x, @SVector [i, j, k]) * y(i, j, k)
# grady((i, j, k)) = ForwardDiff.gradient(_y, @SVector [i, j, k])
# grady((4, 5, 6))
# end

save_vtk(mesh, I)

cm = conservative_metrics(mesh, (5, 14, 5))
i, j, k = 5, 14, 4
cm_i₊½ = conservative_metrics(mesh, (i + 0.5, j, k))
cm_i₋½ = conservative_metrics(mesh, (i - 0.5, j, k))
cm_j₊½ = conservative_metrics(mesh, (i, j + 0.5, k))
cm_j₋½ = conservative_metrics(mesh, (i, j - 0.5, k))
cm_k₊½ = conservative_metrics(mesh, (i, j, k + 0.5))
cm_k₋½ = conservative_metrics(mesh, (i, j, k - 0.5))

cm_i₊½.ξ̂.x - cm_i₋½.ξ̂.x
cm_i₊½.ξ̂.y - cm_i₋½.ξ̂.y
cm_i₊½.ξ̂.z - cm_i₋½.ξ̂.z

cm_k₊½.ζ̂.x - cm_k₋½.ζ̂.x
cm_k₊½.ζ̂.y - cm_k₋½.ζ̂.y
cm_k₊½.ζ̂.z - cm_k₋½.ζ̂.z

I₁ = (cm_i₊½.ξ̂.x - cm_i₋½.ξ̂.x) + (cm_j₊½.η̂.x - cm_j₋½.η̂.x) + (cm_k₊½.ζ̂.x - cm_k₋½.ζ̂.x)
I₂ = (cm_i₊½.ξ̂.y - cm_i₋½.ξ̂.y) + (cm_j₊½.η̂.y - cm_j₋½.η̂.y) + (cm_k₊½.ζ̂.y - cm_k₋½.ζ̂.y)
I₃ = (cm_i₊½.ξ̂.z - cm_i₋½.ξ̂.z) + (cm_j₊½.η̂.z - cm_j₋½.η̂.z) + (cm_k₊½.ζ̂.z - cm_k₋½.ζ̂.z)

I₂ - I₃

i₊½ = metrics(mesh, (5.5, 14, 5), 0)
i₋½ = metrics(mesh, (4.5, 14, 5), 0)

i₊½.ξ.x - i₋½.ξ.x
i₊½.ξ.y - i₋½.ξ.y
i₊½.ξ.z - i₋½.ξ.z

begin
  i, j, k = 5, 14, 4
  m_i₊½ = metrics(mesh, (i + 0.5, j, k), 0)
  m_i₋½ = metrics(mesh, (i - 0.5, j, k), 0)
  m_j₊½ = metrics(mesh, (i, j + 0.5, k), 0)
  m_j₋½ = metrics(mesh, (i, j - 0.5, k), 0)
  m_k₊½ = metrics(mesh, (i, j, k + 0.5), 0)
  m_k₋½ = metrics(mesh, (i, j, k - 0.5), 0)
end