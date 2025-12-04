using Test
using KernelAbstractions
using StructArrays

@testset "SphericalGrid3D Geometry Tests"

begin
  backend = CPU()   # works on all machines

  #-------------------------------------------------------------
  # Simple uniform spherical grid
  #-------------------------------------------------------------
  r = collect(range(1.0, 2.0; length=30))       # 5 radial cells

  θmin = deg2rad(45)
  θmax = deg2rad(105)
  ϕmin = deg2rad(-15)
  ϕmax = deg2rad(15)

  θ = collect(range(θmin, θmax; length=15))       # 4 θ-cells
  ϕ = collect(range(ϕmin, ϕmax; length=25))       # 4 ϕ-cells

  nhalo = 1

  grid = CurvilinearGrids.GridTypes.SphericalGrid3D(r, θ, ϕ, nhalo, backend)

  node = grid.iterators.node
  cell = grid.iterators.cell

  @testset "Interior node xyz correct" begin
    sph = grid.node_coordinates
    cart = grid.cartesian_node_coordinates

    for I in node.domain
      i, j, k = Tuple(I)
      rr = sph.r[i]
      th = sph.θ[j]
      ph = sph.ϕ[k]

      @test isapprox(cart.x[I], rr * sin(th) * cos(ph); atol=1e-12)
      @test isapprox(cart.y[I], rr * sin(th) * sin(ph); atol=1e-12)
      @test isapprox(cart.z[I], rr * cos(th); atol=1e-12)
    end
  end

  @testset "Cell volumes correct" begin
    V = grid.cell_volumes
    sph = grid.node_coordinates

    for I in cell.domain
      i, j, k = Tuple(I)

      r0 = sph.r[i]
      r1 = sph.r[i + 1]
      th0 = sph.θ[j]
      th1 = sph.θ[j + 1]
      ph0 = sph.ϕ[k]
      ph1 = sph.ϕ[k + 1]

      V_true = (1 / 3) * (r1^3 - r0^3) * (cos(th0) - cos(th1)) * (ph1 - ph0)

      @test isapprox(V[I], V_true; rtol=1e-12)
    end
  end

  @testset "Face areas correct" begin
    A = grid.face_areas
    sph = grid.node_coordinates
    rnode, θnode, ϕnode = sph.r, sph.θ, sph.ϕ

    Ai, Aj, Ak = A.i₊½, A.j₊½, A.k₊½

    for I in cell.domain
      i, j, k = Tuple(I)

      Δμ = cos(θnode[j]) - cos(θnode[j + 1])
      Δϕ = ϕnode[k + 1] - ϕnode[k]
      Δr2 = rnode[i + 1]^2 - rnode[i]^2

      Ai_true = rnode[i]^2 * Δμ * Δϕ
      Aj_true = 0.5 * Δr2 * sin(θnode[j]) * Δϕ
      Ak_true = 0.5 * Δr2 * Δμ

      @test isapprox(Ai[I], Ai_true; rtol=1e-12)
      @test isapprox(Aj[I], Aj_true; rtol=1e-12)
      @test isapprox(Ak[I], Ak_true; rtol=1e-12)
    end
  end

  @testset "Centroids correct" begin
    C = grid.centroid_coordinates
    sph = grid.node_coordinates
    rnode, θnode, ϕnode = sph.r, sph.θ, sph.ϕ

    for I in cell.domain
      i, j, k = Tuple(I)

      r₀ = rnode[i]
      r₁ = rnode[i + 1]
      θ₀ = θnode[j]
      θ₁ = θnode[j + 1]
      ϕ₀ = ϕnode[k]
      ϕ₁ = ϕnode[k + 1]

      rc = (3 / 4) * ((r₁^4 - r₀^4) / (r₁^3 - r₀^3))
      θc = acos((cos(θ₀) + cos(θ₁)) / 2)
      ϕc = (ϕ₀ + ϕ₁) / 2

      @test isapprox(C.r[i], rc; rtol=1e-12)

      @test isapprox(C.θ[j], θc; rtol=1e-12)

      @test isapprox(C.ϕ[k], ϕc; rtol=1e-12)
    end
  end

  save_vtk(grid, "spherical_mesh")
end
