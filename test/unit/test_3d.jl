
# @testset "3D Rectangular Mesh Metrics, Conserved Metrics, GCL"
# using CurvilinearGrids
# begin

#   x0, x1 = (0.0, 2.0)
#   y0, y1 = (1.0, 3.0)
#   z0, z1 = (-1.0, 2.0)
#   ni, nj, nk = (40, 80, 120)

#   mesh = rectilinear_grid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :meg6)
#   @code_warntype rectilinear_grid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :meg6)

#   nothing
# end

function rectilinear_nodes_3d(
  x0, x1, y0, y1, z0, z1, ni, nj, nk; halo_coords_included=false
)
  Δx = (x1 - x0) / ni
  Δy = (y1 - y0) / nj
  Δz = (z1 - z0) / nk

  if halo_coords_included
    xnodes = collect(range(x0 + Δx; step=Δx, length=ni + 1))
    ynodes = collect(range(y0 + Δy; step=Δy, length=nj + 1))
    znodes = collect(range(z0 + Δz; step=Δz, length=nk + 1))
  else
    xnodes = collect(range(x0, x1; length=ni + 1))
    ynodes = collect(range(y0, y1; length=nj + 1))
    znodes = collect(range(z0, z1; length=nk + 1))
  end

  x = [xnodes[i] for i in eachindex(xnodes), j in eachindex(ynodes), k in eachindex(znodes)]
  y = [ynodes[j] for i in eachindex(xnodes), j in eachindex(ynodes), k in eachindex(znodes)]
  z = [znodes[k] for i in eachindex(xnodes), j in eachindex(ynodes), k in eachindex(znodes)]
  return x, y, z, Δx, Δy, Δz
end

function assert_rectilinear_metrics_3d(mesh, domain, Δx, Δy, Δz)
  cm = cell_metrics(mesh)
  fm = face_metrics(mesh)

  J = Δx * Δy * Δz
  ξx = inv(Δx)
  ηy = inv(Δy)
  ζz = inv(Δz)
  ξ̂x = J * ξx
  η̂y = J * ηy
  ζ̂z = J * ζz

  @test all(cm.forward[idx].J ≈ J for idx in domain)
  @test all(cm.forward[idx][1, 1] ≈ Δx for idx in domain)
  @test all(cm.forward[idx][2, 2] ≈ Δy for idx in domain)
  @test all(cm.forward[idx][3, 3] ≈ Δz for idx in domain)
  @test all(cm.forward[idx][1, 2] ≈ 0.0 for idx in domain)
  @test all(cm.forward[idx][1, 3] ≈ 0.0 for idx in domain)
  @test all(cm.forward[idx][2, 1] ≈ 0.0 for idx in domain)
  @test all(cm.forward[idx][2, 3] ≈ 0.0 for idx in domain)
  @test all(cm.forward[idx][3, 1] ≈ 0.0 for idx in domain)
  @test all(cm.forward[idx][3, 2] ≈ 0.0 for idx in domain)

  @test all(cm.inverse[idx][1, 1] ≈ ξx for idx in domain)
  @test all(cm.inverse[idx][2, 2] ≈ ηy for idx in domain)
  @test all(cm.inverse[idx][3, 3] ≈ ζz for idx in domain)
  @test all(cm.inverse[idx][1, 2] ≈ 0.0 for idx in domain)
  @test all(cm.inverse[idx][1, 3] ≈ 0.0 for idx in domain)
  @test all(cm.inverse[idx][2, 1] ≈ 0.0 for idx in domain)
  @test all(cm.inverse[idx][2, 3] ≈ 0.0 for idx in domain)
  @test all(cm.inverse[idx][3, 1] ≈ 0.0 for idx in domain)
  @test all(cm.inverse[idx][3, 2] ≈ 0.0 for idx in domain)

  iaxis, jaxis, kaxis = (1, 2, 3)
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)
  k₊½_domain = expand(domain, kaxis, -1)

  @test all(fm[1].forward[idx].J ≈ J for idx in i₊½_domain)
  @test all(fm[2].forward[idx].J ≈ J for idx in j₊½_domain)
  @test all(fm[3].forward[idx].J ≈ J for idx in k₊½_domain)
  @test all(fm[1].forward[idx][1, 1] ≈ Δx for idx in i₊½_domain)
  @test all(fm[2].forward[idx][1, 1] ≈ Δx for idx in j₊½_domain)
  @test all(fm[3].forward[idx][1, 1] ≈ Δx for idx in k₊½_domain)
  @test all(fm[1].forward[idx][2, 2] ≈ Δy for idx in i₊½_domain)
  @test all(fm[2].forward[idx][2, 2] ≈ Δy for idx in j₊½_domain)
  @test all(fm[3].forward[idx][2, 2] ≈ Δy for idx in k₊½_domain)
  @test all(fm[1].forward[idx][3, 3] ≈ Δz for idx in i₊½_domain)
  @test all(fm[2].forward[idx][3, 3] ≈ Δz for idx in j₊½_domain)
  @test all(fm[3].forward[idx][3, 3] ≈ Δz for idx in k₊½_domain)

  @test all(fm[1].inverse[idx][1, 1] ≈ ξx for idx in i₊½_domain)
  @test all(fm[2].inverse[idx][1, 1] ≈ ξx for idx in j₊½_domain)
  @test all(fm[3].inverse[idx][1, 1] ≈ ξx for idx in k₊½_domain)
  @test all(fm[1].inverse[idx][2, 2] ≈ ηy for idx in i₊½_domain)
  @test all(fm[2].inverse[idx][2, 2] ≈ ηy for idx in j₊½_domain)
  @test all(fm[3].inverse[idx][2, 2] ≈ ηy for idx in k₊½_domain)
  @test all(fm[1].inverse[idx][3, 3] ≈ ζz for idx in i₊½_domain)
  @test all(fm[2].inverse[idx][3, 3] ≈ ζz for idx in j₊½_domain)
  @test all(fm[3].inverse[idx][3, 3] ≈ ζz for idx in k₊½_domain)

  @test all(fm[1].conserved[idx][1, 1] ≈ ξ̂x for idx in i₊½_domain)
  @test all(fm[2].conserved[idx][1, 1] ≈ ξ̂x for idx in j₊½_domain)
  @test all(fm[3].conserved[idx][1, 1] ≈ ξ̂x for idx in k₊½_domain)
  @test all(fm[1].conserved[idx][2, 2] ≈ η̂y for idx in i₊½_domain)
  @test all(fm[2].conserved[idx][2, 2] ≈ η̂y for idx in j₊½_domain)
  @test all(fm[3].conserved[idx][2, 2] ≈ η̂y for idx in k₊½_domain)
  @test all(fm[1].conserved[idx][3, 3] ≈ ζ̂z for idx in i₊½_domain)
  @test all(fm[2].conserved[idx][3, 3] ≈ ζ̂z for idx in j₊½_domain)
  @test all(fm[3].conserved[idx][3, 3] ≈ ζ̂z for idx in k₊½_domain)
end

@testset "3D Rectangular Mesh Metrics, Conserved Metrics, GCL" begin
  x0, x1 = (0.0, 2.0)
  y0, y1 = (1, 3)
  z0, z1 = (-1, 2)
  ni, nj, nk = (40, 80, 120)

  x, y, z, Δx, Δy, Δz = rectilinear_nodes_3d(x0, x1, y0, y1, z0, z1, ni, nj, nk)
  mesh = DiscreteGrid(x, y, z, 5)
  domain = mesh.iterators.cell.domain

  assert_rectilinear_metrics_3d(mesh, domain, Δx, Δy, Δz)

  I₁, I₂, I₃ = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )
  domain = mesh.iterators.cell.domain
  gcl_identities = (
    all(abs.(I₁[domain]) .< eps()),
    all(abs.(I₂[domain]) .< eps()),
    all(abs.(I₃[domain]) .< eps()),
  )
  max_vals = (maximum(abs, I₁[domain]), maximum(abs, I₂[domain]), maximum(abs, I₃[domain]))
  @test all(gcl_identities)
end

@testset "Legacy 3D jacobian_matrix matches forward metrics" begin
  x, y, z, _, _, _ = rectilinear_nodes_3d(0.0, 2.0, 1.0, 3.0, -1.0, 2.0, 6, 7, 8)
  mesh = CurvilinearGrid3D(x, y, z, :meg6)
  idx = first(mesh.iterators.cell.domain)
  expected = @SMatrix [
    mesh.cell_center_metrics.x₁.ξ[idx] mesh.cell_center_metrics.x₁.η[idx] mesh.cell_center_metrics.x₁.ζ[idx]
    mesh.cell_center_metrics.x₂.ξ[idx] mesh.cell_center_metrics.x₂.η[idx] mesh.cell_center_metrics.x₂.ζ[idx]
    mesh.cell_center_metrics.x₃.ξ[idx] mesh.cell_center_metrics.x₃.η[idx] mesh.cell_center_metrics.x₃.ζ[idx]
  ]
  @test jacobian_matrix(mesh, idx.I) ≈ expected
end

@testset "3D Wavy Mesh GCL" begin
  using CurvilinearGrids
  using WriteVTK

  function wavy_grid(ni, nj, nk)
    Lx = Ly = Lz = 4.0

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    x = zeros(ni, nj, nk)
    y = zeros(ni, nj, nk)
    z = zeros(ni, nj, nk)
    for k in 1:nk
      for j in 1:nj
        for i in 1:ni
          x[i, j, k] = xmin + Δx0 * ((i - 1) + sinpi((j - 1) * Δy0) * sinpi((k - 1) * Δz0))
          y[i, j, k] = ymin + Δy0 * ((j - 1) + sinpi((k - 1) * Δz0) * sinpi((i - 1) * Δx0))
          z[i, j, k] = zmin + Δz0 * ((k - 1) + sinpi((i - 1) * Δx0) * sinpi((j - 1) * Δy0))
        end
      end
    end

    return (x, y, z)
  end

  ni = nj = nk = 20
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = DiscreteGrid(x, y, z, 5)

  save_vtk(coords(mesh), "wavy3d")

  I₁, I₂, I₃ = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )
  domain = mesh.iterators.cell.domain
  gcl_identities = (
    all(abs.(I₁[domain]) .< 5e-13),
    all(abs.(I₂[domain]) .< 5e-13),
    all(abs.(I₃[domain]) .< 5e-13),
  )
  max_vals = (maximum(abs, I₁[domain]), maximum(abs, I₂[domain]), maximum(abs, I₃[domain]))
  @test all(gcl_identities)
end

@testset "3D Wavy Mesh GCL -- Halo Defined Geometry" begin
  using CurvilinearGrids
  using WriteVTK

  function wavy_grid(ni, nj, nk)
    Lx = Ly = Lz = 4.0

    xmin = -Lx / 2
    ymin = -Ly / 2
    zmin = -Lz / 2

    Δx0 = Lx / ni
    Δy0 = Ly / nj
    Δz0 = Lz / nk

    x = zeros(ni, nj, nk)
    y = zeros(ni, nj, nk)
    z = zeros(ni, nj, nk)
    for k in 1:nk
      for j in 1:nj
        for i in 1:ni
          x[i, j, k] = xmin + Δx0 * ((i - 1) + sinpi((j - 1) * Δy0) * sinpi((k - 1) * Δz0))
          y[i, j, k] = ymin + Δy0 * ((j - 1) + sinpi((k - 1) * Δz0) * sinpi((i - 1) * Δx0))
          z[i, j, k] = zmin + Δz0 * ((k - 1) + sinpi((i - 1) * Δx0) * sinpi((j - 1) * Δy0))
        end
      end
    end

    return (x, y, z)
  end

  ni = nj = nk = 20
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = DiscreteGrid(x, y, z, 5; halo_coords_included=true)

  I₁, I₂, I₃ = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )
  domain = mesh.iterators.cell.domain
  gcl_identities = (
    all(abs.(I₁[domain]) .< eps()),
    all(abs.(I₂[domain]) .< eps()),
    all(abs.(I₃[domain]) .< eps()),
  )
  max_vals = (maximum(abs, I₁[domain]), maximum(abs, I₂[domain]), maximum(abs, I₃[domain]))
  @test all(gcl_identities)
end

@testset "3D Sphere Sector, Symmetric Conservative Metrics" begin
  r0, r1 = (1, 3)
  (θ0, θ1) = deg2rad.((35, 180 - 35))
  (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))

  ni, nj, nk = (20, 20, 20)
  params = (; r0, θ0, ϕ0, Δr=(r1 - r0) / ni, Δθ=(θ1 - θ0) / nj, Δϕ=(ϕ1 - ϕ0) / nk)

  x(t, ξ, η, ζ, p) =
    (p.r0 + (ξ - 1) * p.Δr) * sin(p.θ0 + (η - 1) * p.Δθ) * cos(p.ϕ0 + (ζ - 1) * p.Δϕ)
  y(t, ξ, η, ζ, p) =
    (p.r0 + (ξ - 1) * p.Δr) * sin(p.θ0 + (η - 1) * p.Δθ) * sin(p.ϕ0 + (ζ - 1) * p.Δϕ)
  z(t, ξ, η, ζ, p) = (p.r0 + (ξ - 1) * p.Δr) * cos(p.θ0 + (η - 1) * p.Δθ)

  mesh = MappedGrid(x, y, z, params, (ni, nj, nk), 5)
  save_vtk(coords(mesh), "sphere_sector_3d")

  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  I₁, I₂, I₃ = CurvilinearGrids.GridTypes.gcl(
    face_metrics(mesh), mesh.iterators.cell.domain
  )
  domain = mesh.iterators.cell.domain
  gcl_identities = (
    all(abs.(I₁[domain]) .< 5e-13),
    all(abs.(I₂[domain]) .< 5e-13),
    all(abs.(I₃[domain]) .< 5e-13),
  )
  max_vals = (maximum(abs, I₁[domain]), maximum(abs, I₂[domain]), maximum(abs, I₃[domain]))
  # @show gcl_identities, max_vals
  @test all(gcl_identities)

  nothing
end
