using Test
using KernelAbstractions
using StructArrays
using CurvilinearGrids
using StaticArrays

# @testset "SphericalGrid3D Geometry Tests" begin
#   backend = CPU()   # works on all machines

#   #-------------------------------------------------------------
#   # Simple uniform spherical grid
#   #-------------------------------------------------------------
#   r = collect(range(1.0, 2.0; length=30))       # 5 radial cells

#   θmin = deg2rad(45)
#   θmax = deg2rad(105)
#   ϕmin = deg2rad(-15)
#   ϕmax = deg2rad(15)

#   θ = collect(range(θmin, θmax; length=15))       # 4 θ-cells
#   ϕ = collect(range(ϕmin, ϕmax; length=25))       # 4 ϕ-cells

#   nhalo = 1

#   grid = CurvilinearGrids.GridTypes.SphericalGrid3D(r, θ, ϕ, nhalo, backend)

#   node = grid.iterators.node
#   cell = grid.iterators.cell

#   @testset "Interior node xyz correct" begin
#     sph = grid.node_coordinates
#     cart = grid.cartesian_node_coordinates

#     for I in node.domain
#       i, j, k = Tuple(I)
#       rr = sph.r[i]
#       th = sph.θ[j]
#       ph = sph.ϕ[k]

#       @test isapprox(cart.x[I], rr * sin(th) * cos(ph); atol=1e-12)
#       @test isapprox(cart.y[I], rr * sin(th) * sin(ph); atol=1e-12)
#       @test isapprox(cart.z[I], rr * cos(th); atol=1e-12)
#     end
#   end

#   @testset "Cell volumes correct" begin
#     V = grid.cell_volumes
#     sph = grid.node_coordinates

#     for I in cell.domain
#       i, j, k = Tuple(I)

#       r0 = sph.r[i]
#       r1 = sph.r[i + 1]
#       th0 = sph.θ[j]
#       th1 = sph.θ[j + 1]
#       ph0 = sph.ϕ[k]
#       ph1 = sph.ϕ[k + 1]

#       V_true = (1 / 3) * (r1^3 - r0^3) * (cos(th0) - cos(th1)) * (ph1 - ph0)

#       @test isapprox(V[I], V_true; rtol=1e-12)
#     end
#   end

#   @testset "Face areas correct" begin
#     A = grid.face_areas
#     sph = grid.node_coordinates
#     rnode, θnode, ϕnode = sph.r, sph.θ, sph.ϕ

#     Ai, Aj, Ak = A.i₊½, A.j₊½, A.k₊½

#     for I in cell.domain
#       i, j, k = Tuple(I)

#       Δμ = cos(θnode[j]) - cos(θnode[j + 1])
#       Δϕ = ϕnode[k + 1] - ϕnode[k]
#       Δr2 = rnode[i + 1]^2 - rnode[i]^2

#       Ai_true = rnode[i]^2 * Δμ * Δϕ
#       Aj_true = 0.5 * Δr2 * sin(θnode[j]) * Δϕ
#       Ak_true = 0.5 * Δr2 * Δμ

#       @test isapprox(Ai[I], Ai_true; rtol=1e-12)
#       @test isapprox(Aj[I], Aj_true; rtol=1e-12)
#       @test isapprox(Ak[I], Ak_true; rtol=1e-12)
#     end
#   end

#   @testset "Centroids correct" begin
#     C = grid.centroid_coordinates
#     sph = grid.node_coordinates
#     rnode, θnode, ϕnode = sph.r, sph.θ, sph.ϕ

#     for I in cell.domain
#       i, j, k = Tuple(I)

#       r₀ = rnode[i]
#       r₁ = rnode[i + 1]
#       θ₀ = θnode[j]
#       θ₁ = θnode[j + 1]
#       ϕ₀ = ϕnode[k]
#       ϕ₁ = ϕnode[k + 1]

#       rc = (3 / 4) * ((r₁^4 - r₀^4) / (r₁^3 - r₀^3))
#       θc = acos((cos(θ₀) + cos(θ₁)) / 2)
#       ϕc = (ϕ₀ + ϕ₁) / 2

#       @test isapprox(C.r[i], rc; rtol=1e-12)

#       @test isapprox(C.θ[j], θc; rtol=1e-12)

#       @test isapprox(C.ϕ[k], ϕc; rtol=1e-12)
#     end
#   end

#   save_vtk(grid, "spherical_mesh")
# end

# analytic scalar field:
#   T(r,θ,ϕ) = r^2 + 2cosθ + 3sinϕ
# gradient:
#   ∂T/∂r = 2r
#   ∂T/∂θ = -2sinθ
#   ∂T/∂ϕ = 3cosϕ
#
# analytic vector field:
#   A = (A_r, A_θ, A_ϕ) = (r, 0, 0)
# divergence:
#   ∇·A = 3
# curl:
#   ∇×A = 0

@testset "SphericalGrid cell-center operators" begin
  backend = KernelAbstractions.CPU()
  nhalo = 2

  θmid = π / 2
  ϕmid = 0.0
  nr, ntheta, nphi = (50, 50, 50)

  r_nodes = range(1.0, 2.0; length=nr) |> collect
  θ_nodes = range(θmid - 0.1, θmid + 0.1; length=ntheta) |> collect
  ϕ_nodes = range(ϕmid - 0.1, ϕmid + 0.1; length=nphi) |> collect

  grid = SphericalGrid3D(
    r_nodes, θ_nodes, ϕ_nodes, nhalo, backend; halo_coords_included=true
  )

  # fill analytic fields on interior domain
  cell_dom = grid.iterators.cell.domain
  full_dom = grid.iterators.cell.full

  V = grid.cell_volumes
  rc = grid.centroid_coordinates.r
  θc = grid.centroid_coordinates.θ
  ϕc = grid.centroid_coordinates.ϕ

  # allocate cell-centered scalar T and vector A
  T = similar(V)
  dT_dr = similar(V)
  dT_dθ = similar(V)
  dT_dϕ = similar(V)
  dT_dr_analytic = similar(V)
  dT_dθ_analytic = similar(V)
  dT_dϕ_analytic = similar(V)
  A_r = similar(V)
  A_θ = similar(V)
  A_ϕ = similar(V)

  for I in full_dom
    i, j, k = I.I
    r = rc[i]
    θ = θc[j]
    ϕ = ϕc[k]

    T[I] = r^2 + 2 * cos(θ) + 3 * sin(ϕ)
    A_r[I] = r
    A_θ[I] = 0.0
    A_ϕ[I] = 0.0
  end

  # -------- gradient test --------
  max_grad_err = 0.0
  for I in cell_dom
    i, j, k = I.I
    r = rc[i]
    θ = θc[j]
    ϕ = ϕc[k]

    g_num = cell_center_gradient(grid, T, I)

    dTdr, dTdθ, dTdϕ = g_num
    dT_dr[I] = dTdr
    dT_dθ[I] = dTdθ
    dT_dϕ[I] = dTdϕ

    g_true = @SVector [2r, -2 * sin(θ) / r, 3 * cos(ϕ) / (r * sin(θ))]

    dTdr_true, dTdθ_true, dTdϕ_true = g_true
    dT_dr_analytic[I] = dTdr_true
    dT_dθ_analytic[I] = dTdθ_true
    dT_dϕ_analytic[I] = dTdϕ_true

    err = maximum(abs.(g_num - g_true))
    max_grad_err = max(max_grad_err, err)
  end
  @info "Max gradient error ≈ $max_grad_err"
  @test max_grad_err < 0.015

  # -------- divergence test --------
  max_div_err = 0.0
  for I in cell_dom
    div_num = cell_center_divergence(grid, (A_r, A_θ, A_ϕ), I)
    div_true = 3.0
    max_div_err = max(max_div_err, abs(div_num - div_true))
  end
  @info "Max divergence error ≈ $max_div_err"
  @test max_div_err < 1e-2

  # -------- curl test (should be ~0) --------
  max_curl_norm = 0.0
  for I in cell_dom
    curl_num = cell_center_curl(grid, (A_r, A_θ, A_ϕ), I)
    max_curl_norm = max(max_curl_norm, maximum(abs.(curl_num)))
  end
  @info "Max curl error ≈ $max_curl_norm"
  @test max_curl_norm < 1e-1

  CurvilinearGrids.save_vtk(
    grid,
    "sphere";
    extra_cell_data=(;
      T, A_r, A_θ, A_ϕ, dT_dr, dT_dθ, dT_dϕ, dT_dr_analytic, dT_dθ_analytic, dT_dϕ_analytic
    ), # add data to plot
  )
end
