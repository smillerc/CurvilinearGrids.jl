using Test

using CurvilinearGrids
using StaticArrays
using Test
using BenchmarkTools
using WriteVTK

function save_vtk(mesh)
  fn = "wavy"
  @info "Writing to $fn.vti"

  xyz_n = CurvilinearGrids.coords(mesh)
  domain = mesh.iterators.cell.domain

  @views begin
    J = [m.J for m in mesh.cell_center_metrics[domain]]
    ξx = [m.ξx for m in mesh.cell_center_metrics[domain]]
    ξ̂xᵢ₊½ = [m.ξ̂x for m in mesh.edge_metrics.i₊½[domain]]
    ηy = [m.ηy for m in mesh.cell_center_metrics[domain]]
  end

  vtk_grid(fn, xyz_n) do vtk

    # w1 = [x[1] for idx in eachindex(f.char_state.W[ilo:ihi, jlo:jhi]

    vtk["J"] = J
    vtk["ξx"] = ξx
    vtk["ξ̂xᵢ₊½"] = ξ̂xᵢ₊½
    vtk["ηy"] = ηy
  end
end

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
mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo);

save_vtk(mesh);

ξ̂xᵢ₊½ = [m.ξ̂x for m in mesh.edge_metrics.k₊½]
size(ξ̂xᵢ₊½)
using Plots
heatmap(ξ̂xᵢ₊½[:, :, 3])