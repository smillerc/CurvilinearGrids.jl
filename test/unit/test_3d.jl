
@testset "3D Rectangular Mesh Metrics, Conserved Metrics, GCL" begin
  include("common.jl")

  x0, x1 = (0.0, 2.0)
  y0, y1 = (1, 3)
  z0, z1 = (-1, 2)
  ni, nj, nk = (40, 80, 120)

  mesh = RectlinearGrid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :meg6)
  domain = mesh.iterators.cell.domain

  cell_volume = 0.05 * 0.025 * 0.025

  @test all(mesh.cell_center_metrics.forward.J[domain] .≈ cell_volume)

  @test all(mesh.cell_center_metrics.forward.x₁.ξ[domain] .≈ 0.05)
  @test all(mesh.cell_center_metrics.forward.x₂.ξ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₃.ξ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₁.η[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₂.η[domain] .≈ 0.025)
  @test all(mesh.cell_center_metrics.forward.x₃.η[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₁.ζ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₂.ζ[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.forward.x₃.ζ[domain] .≈ 0.025)

  @test all(mesh.cell_center_metrics.inverse_normalized.ξ̂.x₁[domain] .≈ 0.000625)
  @test all(mesh.cell_center_metrics.inverse_normalized.ξ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse_normalized.ξ̂.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse_normalized.η̂.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse_normalized.η̂.x₂[domain] .≈ 0.00125)
  @test all(mesh.cell_center_metrics.inverse_normalized.η̂.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse_normalized.ζ̂.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse_normalized.ζ̂.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse_normalized.ζ̂.x₃[domain] .≈ 0.00125)

  @test all(mesh.cell_center_metrics.inverse.ξ.x₁[domain] .≈ 20.0)
  @test all(mesh.cell_center_metrics.inverse.ξ.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.ξ.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.η.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.η.x₂[domain] .≈ 40.0)
  @test all(mesh.cell_center_metrics.inverse.η.x₃[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.ζ.x₁[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.ζ.x₂[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.ζ.x₃[domain] .≈ 40.0)
  @test all(mesh.cell_center_metrics.inverse.ξ.t[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.η.t[domain] .≈ 0.0)
  @test all(mesh.cell_center_metrics.inverse.ζ.t[domain] .≈ 0.0)

  iaxis, jaxis, kaxis = (1, 2, 3)
  i₊½_domain = expand(domain, iaxis, -1)
  j₊½_domain = expand(domain, jaxis, -1)
  k₊½_domain = expand(domain, kaxis, -1)

  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₁[i₊½_domain] .≈ 0.000625)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ξ̂.x₁[j₊½_domain] .≈ 0.000625)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ξ̂.x₁[k₊½_domain] .≈ 0.000625)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₂[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ξ̂.x₂[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ξ̂.x₂[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.x₃[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ξ̂.x₃[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ξ̂.x₃[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.η̂.x₁[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.η̂.x₁[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.η̂.x₁[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.η̂.x₂[i₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.η̂.x₂[j₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.η̂.x₂[k₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.η̂.x₃[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.η̂.x₃[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.η̂.x₃[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ζ̂.x₁[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ζ̂.x₁[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ζ̂.x₁[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ζ̂.x₂[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ζ̂.x₂[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ζ̂.x₂[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ζ̂.x₃[i₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ζ̂.x₃[j₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ζ̂.x₃[k₊½_domain] .≈ 0.00125)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ξ̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ξ̂.t[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ξ̂.t[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.η̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.η̂.t[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.η̂.t[k₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.i₊½.ζ̂.t[i₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.j₊½.ζ̂.t[j₊½_domain] .≈ 0.0)
  @test all(mesh.edge_metrics.inverse_normalized.k₊½.ζ̂.t[k₊½_domain] .≈ 0.0)

  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  ϵ = 5e-15
  em = mesh.edge_metrics.inverse_normalized

  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (em.i₊½.ξ̂.x₁[i, j, k] - em.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (em.j₊½.η̂.x₁[i, j, k] - em.j₊½.η̂.x₁[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₁[i, j, k] - em.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂ = (
      (em.i₊½.ξ̂.x₂[i, j, k] - em.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (em.j₊½.η̂.x₂[i, j, k] - em.j₊½.η̂.x₂[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₂[i, j, k] - em.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃ = (
      (em.i₊½.ξ̂.x₃[i, j, k] - em.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (em.j₊½.η̂.x₃[i, j, k] - em.j₊½.η̂.x₃[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₃[i, j, k] - em.k₊½.ζ̂.x₃[i, j, k - 1])
    )

    I₁_passes = abs(I₁) < ϵ
    I₂_passes = abs(I₂) < ϵ
    I₃_passes = abs(I₃) < ϵ
    if !(I₁_passes && I₂_passes && I₃_passes)
      @show I₁ I₂ I₃
      break
    end
  end
  @test I₁_passes
  @test I₂_passes
  @test I₃_passes
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

  ξ = 1
  η = 2
  ζ = 3

  ni = nj = nk = 20
  x, y, z = wavy_grid(ni, nj, nk)
  mesh = CurvilinearGrid3D(x, y, z, :meg6)

  save_vtk(mesh, "wavy3d")

  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  ϵ = 5e-14
  # ϵ = 5e-15
  em = mesh.edge_metrics.inverse_normalized

  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (em.i₊½.ξ̂.x₁[i, j, k] - em.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (em.j₊½.η̂.x₁[i, j, k] - em.j₊½.η̂.x₁[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₁[i, j, k] - em.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂ = (
      (em.i₊½.ξ̂.x₂[i, j, k] - em.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (em.j₊½.η̂.x₂[i, j, k] - em.j₊½.η̂.x₂[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₂[i, j, k] - em.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃ = (
      (em.i₊½.ξ̂.x₃[i, j, k] - em.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (em.j₊½.η̂.x₃[i, j, k] - em.j₊½.η̂.x₃[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₃[i, j, k] - em.k₊½.ζ̂.x₃[i, j, k - 1])
    )

    I₁_passes = abs(I₁) < ϵ
    I₂_passes = abs(I₂) < ϵ
    I₃_passes = abs(I₃) < ϵ
    if !(I₁_passes && I₂_passes && I₃_passes)
      @show I₁ I₂ I₃
      break
    end
  end
  @test I₁_passes
  @test I₂_passes
  @test I₃_passes
end

@testset "3D Sphere Sector, Symmetric Conservative Metrics" begin
  include("common.jl")

  r0, r1 = (1, 3)
  (θ0, θ1) = deg2rad.((35, 180 - 35))
  (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))

  ni, nj, nk = (20, 20, 20)
  mesh = RThetaPhiGrid((r0, θ0, ϕ0), (r1, θ1, ϕ1), (ni, nj, nk), :meg6_symmetric)
  # mesh = RThetaPhiGrid((r0, θ0, ϕ0), (r1, θ1, ϕ1), (ni, nj, nk), :meg6)
  save_vtk(mesh, "sphere_sector_3d")

  I₁_passes = true
  I₂_passes = true
  I₃_passes = true

  ϵ = 5e-13
  # ϵ = 5e-15
  em = mesh.edge_metrics.inverse_normalized

  I₁_max = -Inf
  I₂_max = -Inf
  I₃_max = -Inf

  for idx in mesh.iterators.cell.domain
    i, j, k = idx.I
    I₁ = (
      (em.i₊½.ξ̂.x₁[i, j, k] - em.i₊½.ξ̂.x₁[i - 1, j, k]) +
      (em.j₊½.η̂.x₁[i, j, k] - em.j₊½.η̂.x₁[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₁[i, j, k] - em.k₊½.ζ̂.x₁[i, j, k - 1])
    )
    I₂ = (
      (em.i₊½.ξ̂.x₂[i, j, k] - em.i₊½.ξ̂.x₂[i - 1, j, k]) +
      (em.j₊½.η̂.x₂[i, j, k] - em.j₊½.η̂.x₂[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₂[i, j, k] - em.k₊½.ζ̂.x₂[i, j, k - 1])
    )
    I₃ = (
      (em.i₊½.ξ̂.x₃[i, j, k] - em.i₊½.ξ̂.x₃[i - 1, j, k]) +
      (em.j₊½.η̂.x₃[i, j, k] - em.j₊½.η̂.x₃[i, j - 1, k]) +
      (em.k₊½.ζ̂.x₃[i, j, k] - em.k₊½.ζ̂.x₃[i, j, k - 1])
    )

    I₁_max = max(I₁_max, abs(I₁))
    I₂_max = max(I₂_max, abs(I₂))
    I₃_max = max(I₃_max, abs(I₃))

    I₁_passes = abs(I₁) < ϵ
    I₂_passes = abs(I₂) < ϵ
    I₃_passes = abs(I₃) < ϵ
    if !(I₁_passes && I₂_passes && I₃_passes)
      @show I₁ I₂ I₃
      break
    end
  end

  @test I₁_passes
  @test I₂_passes
  @test I₃_passes

  @test all(
    abs.(mesh.cell_center_metrics.inverse.ζ.x₃[mesh.iterators.cell.domain]) .< 5e-12
  )

  # @show I₁_max
  # @show I₂_max
  # @show I₃_max

  nothing
end
