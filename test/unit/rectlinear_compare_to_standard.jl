using Test, Random, CurvilinearGrids
using CurvilinearGrids.RectilinearArrays

@testset "Rectilinear3D Equivalence Test" begin
  x = (1:0.5:5) .^ 2
  y = (0:0.3:6) .^ 2
  z = (2:0.2:7) .^ 2
  mesh_rect = RectilinearGrid3D(x, y, z, :MEG6)
  mesh_curv = rectilinear_grid(x, y, z, :MEG6)

  # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.

  domain = mesh_rect.iterators.cell.domain

  # Cell center metrics
  @test isapprox(
    mesh_rect.cell_center_metrics.J[domain],
    mesh_curv.cell_center_metrics.J[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.η[domain],
    mesh_curv.cell_center_metrics.x₁.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ξ[domain],
    mesh_curv.cell_center_metrics.x₁.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ζ[domain],
    mesh_curv.cell_center_metrics.x₁.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.η[domain],
    mesh_curv.cell_center_metrics.x₂.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ξ[domain],
    mesh_curv.cell_center_metrics.x₂.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ζ[domain],
    mesh_curv.cell_center_metrics.x₂.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.η[domain],
    mesh_curv.cell_center_metrics.x₃.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.ξ[domain],
    mesh_curv.cell_center_metrics.x₃.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.ζ[domain],
    mesh_curv.cell_center_metrics.x₃.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η.t[domain],
    mesh_curv.cell_center_metrics.η.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₁[domain],
    mesh_curv.cell_center_metrics.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₂[domain],
    mesh_curv.cell_center_metrics.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₃[domain],
    mesh_curv.cell_center_metrics.η.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.t[domain],
    mesh_curv.cell_center_metrics.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₁[domain],
    mesh_curv.cell_center_metrics.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₂[domain],
    mesh_curv.cell_center_metrics.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₃[domain],
    mesh_curv.cell_center_metrics.η̂.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.t[domain],
    mesh_curv.cell_center_metrics.ξ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₁[domain],
    mesh_curv.cell_center_metrics.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₂[domain],
    mesh_curv.cell_center_metrics.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₃[domain],
    mesh_curv.cell_center_metrics.ξ.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.t[domain],
    mesh_curv.cell_center_metrics.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₃[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.t[domain],
    mesh_curv.cell_center_metrics.ζ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₁[domain],
    mesh_curv.cell_center_metrics.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₂[domain],
    mesh_curv.cell_center_metrics.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₃[domain],
    mesh_curv.cell_center_metrics.ζ.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.t[domain],
    mesh_curv.cell_center_metrics.ζ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₃[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₃[domain],
    atol=1e-10,
  )

  # Edge metrics
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₃[domain],
    mesh_curv.edge_metrics.i₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.t[domain],
    mesh_curv.edge_metrics.i₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₃[domain],
    mesh_curv.edge_metrics.j₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.t[domain],
    mesh_curv.edge_metrics.j₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₁[domain],
    mesh_curv.edge_metrics.k₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₂[domain],
    mesh_curv.edge_metrics.k₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₃[domain],
    mesh_curv.edge_metrics.k₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.t[domain],
    mesh_curv.edge_metrics.k₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.t[domain],
    atol=1e-10,
  )
end

@testset "Rectilinear3D Equivalence Test - Semi-Uniform" begin
  (x0, y0, z0) = (0, 1, 2)
  (x1, y1, z1) = (12, 13, 14)
  (ni, nj, nk) = (20, 30, 40)

  mesh_rect = RectilinearGrid3D((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :MEG6)
  mesh_curv = rectilinear_grid((x0, y0, z0), (x1, y1, z1), (ni, nj, nk), :MEG6)

  # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.

  domain = mesh_rect.iterators.cell.domain

  # Cell center metrics
  @test isapprox(
    mesh_rect.cell_center_metrics.J[domain],
    mesh_curv.cell_center_metrics.J[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.η[domain],
    mesh_curv.cell_center_metrics.x₁.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ξ[domain],
    mesh_curv.cell_center_metrics.x₁.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ζ[domain],
    mesh_curv.cell_center_metrics.x₁.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.η[domain],
    mesh_curv.cell_center_metrics.x₂.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ξ[domain],
    mesh_curv.cell_center_metrics.x₂.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ζ[domain],
    mesh_curv.cell_center_metrics.x₂.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.η[domain],
    mesh_curv.cell_center_metrics.x₃.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.ξ[domain],
    mesh_curv.cell_center_metrics.x₃.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.ζ[domain],
    mesh_curv.cell_center_metrics.x₃.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η.t[domain],
    mesh_curv.cell_center_metrics.η.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₁[domain],
    mesh_curv.cell_center_metrics.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₂[domain],
    mesh_curv.cell_center_metrics.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₃[domain],
    mesh_curv.cell_center_metrics.η.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.t[domain],
    mesh_curv.cell_center_metrics.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₁[domain],
    mesh_curv.cell_center_metrics.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₂[domain],
    mesh_curv.cell_center_metrics.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₃[domain],
    mesh_curv.cell_center_metrics.η̂.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.t[domain],
    mesh_curv.cell_center_metrics.ξ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₁[domain],
    mesh_curv.cell_center_metrics.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₂[domain],
    mesh_curv.cell_center_metrics.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₃[domain],
    mesh_curv.cell_center_metrics.ξ.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.t[domain],
    mesh_curv.cell_center_metrics.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₃[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.t[domain],
    mesh_curv.cell_center_metrics.ζ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₁[domain],
    mesh_curv.cell_center_metrics.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₂[domain],
    mesh_curv.cell_center_metrics.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₃[domain],
    mesh_curv.cell_center_metrics.ζ.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.t[domain],
    mesh_curv.cell_center_metrics.ζ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₃[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₃[domain],
    atol=1e-10,
  )

  # Edge metrics
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₃[domain],
    mesh_curv.edge_metrics.i₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.t[domain],
    mesh_curv.edge_metrics.i₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₃[domain],
    mesh_curv.edge_metrics.j₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.t[domain],
    mesh_curv.edge_metrics.j₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₁[domain],
    mesh_curv.edge_metrics.k₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₂[domain],
    mesh_curv.edge_metrics.k₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₃[domain],
    mesh_curv.edge_metrics.k₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.t[domain],
    mesh_curv.edge_metrics.k₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.t[domain],
    atol=1e-10,
  )
end

@testset "Uniform3D Equivalence Test" begin
  (x0, y0, z0) = (0, 1, 2)
  (x1, y1, z1) = (12, 13, 14)
  ∂x = 0.5
  x = x0:∂x:x1
  y = y0:∂x:y1
  z = z0:∂x:z1

  mesh_rect = UniformGrid3D((x0, y0, z0), (x1, y1, z1), ∂x, :MEG6)
  mesh_curv = rectilinear_grid(x, y, z, :MEG6)

  # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.

  domain = mesh_rect.iterators.cell.domain

  # Cell center metrics
  @test isapprox(
    mesh_rect.cell_center_metrics.J[domain],
    mesh_curv.cell_center_metrics.J[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.η[domain],
    mesh_curv.cell_center_metrics.x₁.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ξ[domain],
    mesh_curv.cell_center_metrics.x₁.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ζ[domain],
    mesh_curv.cell_center_metrics.x₁.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.η[domain],
    mesh_curv.cell_center_metrics.x₂.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ξ[domain],
    mesh_curv.cell_center_metrics.x₂.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ζ[domain],
    mesh_curv.cell_center_metrics.x₂.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.η[domain],
    mesh_curv.cell_center_metrics.x₃.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.ξ[domain],
    mesh_curv.cell_center_metrics.x₃.ξ[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₃.ζ[domain],
    mesh_curv.cell_center_metrics.x₃.ζ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η.t[domain],
    mesh_curv.cell_center_metrics.η.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₁[domain],
    mesh_curv.cell_center_metrics.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₂[domain],
    mesh_curv.cell_center_metrics.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₃[domain],
    mesh_curv.cell_center_metrics.η.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.t[domain],
    mesh_curv.cell_center_metrics.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₁[domain],
    mesh_curv.cell_center_metrics.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₂[domain],
    mesh_curv.cell_center_metrics.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₃[domain],
    mesh_curv.cell_center_metrics.η̂.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.t[domain],
    mesh_curv.cell_center_metrics.ξ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₁[domain],
    mesh_curv.cell_center_metrics.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₂[domain],
    mesh_curv.cell_center_metrics.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₃[domain],
    mesh_curv.cell_center_metrics.ξ.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.t[domain],
    mesh_curv.cell_center_metrics.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₃[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.t[domain],
    mesh_curv.cell_center_metrics.ζ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₁[domain],
    mesh_curv.cell_center_metrics.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₂[domain],
    mesh_curv.cell_center_metrics.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ.x₃[domain],
    mesh_curv.cell_center_metrics.ζ.x₃[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.t[domain],
    mesh_curv.cell_center_metrics.ζ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ζ̂.x₃[domain],
    mesh_curv.cell_center_metrics.ζ̂.x₃[domain],
    atol=1e-10,
  )

  # Edge metrics
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₃[domain],
    mesh_curv.edge_metrics.i₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.t[domain],
    mesh_curv.edge_metrics.i₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ζ̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₃[domain],
    mesh_curv.edge_metrics.j₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.t[domain],
    mesh_curv.edge_metrics.j₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ζ̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ξ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₁[domain],
    mesh_curv.edge_metrics.k₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₂[domain],
    mesh_curv.edge_metrics.k₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η.x₃[domain],
    mesh_curv.edge_metrics.k₊½.η.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ζ.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.k₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.η̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.η̂.t[domain],
    mesh_curv.edge_metrics.k₊½.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₁[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₂[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.x₃[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.x₃[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.k₊½.ζ̂.t[domain],
    mesh_curv.edge_metrics.k₊½.ζ̂.t[domain],
    atol=1e-10,
  )
end

@testset "Rectilinear2D Equivalence Test" begin
  Random.seed!(1236)

  x = cumsum(rand(20))
  y = cumsum(rand(30))
  mesh_rect = RectilinearGrid2D(x, y, :MEG6)
  mesh_curv = rectilinear_grid(x, y, :MEG6)

  # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.

  domain = mesh_rect.iterators.cell.domain

  # Cell center metrics
  @test isapprox(
    mesh_rect.cell_center_metrics.J[domain],
    mesh_curv.cell_center_metrics.J[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.η[domain],
    mesh_curv.cell_center_metrics.x₁.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ξ[domain],
    mesh_curv.cell_center_metrics.x₁.ξ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.η[domain],
    mesh_curv.cell_center_metrics.x₂.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ξ[domain],
    mesh_curv.cell_center_metrics.x₂.ξ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η.t[domain],
    mesh_curv.cell_center_metrics.η.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₁[domain],
    mesh_curv.cell_center_metrics.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₂[domain],
    mesh_curv.cell_center_metrics.η.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.t[domain],
    mesh_curv.cell_center_metrics.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₁[domain],
    mesh_curv.cell_center_metrics.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₂[domain],
    mesh_curv.cell_center_metrics.η̂.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.t[domain],
    mesh_curv.cell_center_metrics.ξ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₁[domain],
    mesh_curv.cell_center_metrics.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₂[domain],
    mesh_curv.cell_center_metrics.ξ.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.t[domain],
    mesh_curv.cell_center_metrics.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₂[domain],
    atol=1e-10,
  )

  # Edge metrics
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.t[domain],
    mesh_curv.edge_metrics.i₊½.η̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.t[domain],
    mesh_curv.edge_metrics.j₊½.η̂.t[domain],
    atol=1e-10,
  )
end

@testset "Rectilinear2D Equivalence Test - Semi-Uniform" begin
  Random.seed!(1236)

  (x0, y0) = (0, 2)
  (x1, y1) = (13, 14)
  (ni, nj) = (40, 52)

  mesh_rect = RectilinearGrid2D((x0, y0), (x1, y1), (ni, nj), :MEG6)
  mesh_curv = rectilinear_grid((x0, y0), (x1, y1), (ni, nj), :MEG6)

  # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.

  domain = mesh_rect.iterators.cell.domain

  # Cell center metrics
  @test isapprox(
    mesh_rect.cell_center_metrics.J[domain],
    mesh_curv.cell_center_metrics.J[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.η[domain],
    mesh_curv.cell_center_metrics.x₁.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₁.ξ[domain],
    mesh_curv.cell_center_metrics.x₁.ξ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.η[domain],
    mesh_curv.cell_center_metrics.x₂.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.x₂.ξ[domain],
    mesh_curv.cell_center_metrics.x₂.ξ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η.t[domain],
    mesh_curv.cell_center_metrics.η.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₁[domain],
    mesh_curv.cell_center_metrics.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η.x₂[domain],
    mesh_curv.cell_center_metrics.η.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.t[domain],
    mesh_curv.cell_center_metrics.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₁[domain],
    mesh_curv.cell_center_metrics.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.η̂.x₂[domain],
    mesh_curv.cell_center_metrics.η̂.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.t[domain],
    mesh_curv.cell_center_metrics.ξ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₁[domain],
    mesh_curv.cell_center_metrics.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ.x₂[domain],
    mesh_curv.cell_center_metrics.ξ.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.t[domain],
    mesh_curv.cell_center_metrics.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.cell_center_metrics.ξ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₂[domain],
    atol=1e-10,
  )

  # Edge metrics
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.i₊½.η̂.t[domain],
    mesh_curv.edge_metrics.i₊½.η̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_rect.edge_metrics.j₊½.η̂.t[domain],
    mesh_curv.edge_metrics.j₊½.η̂.t[domain],
    atol=1e-10,
  )
end

@testset "Uniform2D Equivalence Test" begin
  (x0, y0) = (0, 2)
  (x1, y1) = (14, 15)
  ∂x = 0.5
  x = Vector(x0:∂x:x1)
  y = Vector(y0:∂x:y1)
  mesh_uni = UniformGrid2D((x0, y0), (x1, y1), ∂x, :MEG6)
  mesh_curv = rectilinear_grid(x, y, :MEG6)

  # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.

  domain = mesh_uni.iterators.cell.domain

  # Cell center metrics
  @test isapprox(
    mesh_uni.cell_center_metrics.J[domain],
    mesh_curv.cell_center_metrics.J[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.x₁.η[domain],
    mesh_curv.cell_center_metrics.x₁.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.x₁.ξ[domain],
    mesh_curv.cell_center_metrics.x₁.ξ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.x₂.η[domain],
    mesh_curv.cell_center_metrics.x₂.η[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.x₂.ξ[domain],
    mesh_curv.cell_center_metrics.x₂.ξ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.η.t[domain],
    mesh_curv.cell_center_metrics.η.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.η.x₁[domain],
    mesh_curv.cell_center_metrics.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.η.x₂[domain],
    mesh_curv.cell_center_metrics.η.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.η̂.t[domain],
    mesh_curv.cell_center_metrics.η̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.η̂.x₁[domain],
    mesh_curv.cell_center_metrics.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.η̂.x₂[domain],
    mesh_curv.cell_center_metrics.η̂.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.ξ.t[domain],
    mesh_curv.cell_center_metrics.ξ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.ξ.x₁[domain],
    mesh_curv.cell_center_metrics.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.ξ.x₂[domain],
    mesh_curv.cell_center_metrics.ξ.x₂[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.ξ̂.t[domain],
    mesh_curv.cell_center_metrics.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.ξ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.ξ̂.x₂[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₂[domain],
    atol=1e-10,
  )

  # Edge metrics
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.η.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.η.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.i₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.η̂.t[domain],
    mesh_curv.edge_metrics.i₊½.η̂.t[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.edge_metrics.j₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.ξ.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.η.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.η.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.ξ̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.j₊½.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.η̂.x₁[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.η̂.x₂[domain],
    mesh_curv.edge_metrics.j₊½.η̂.x₂[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.j₊½.η̂.t[domain],
    mesh_curv.edge_metrics.j₊½.η̂.t[domain],
    atol=1e-10,
  )
end

@testset "Uniform1D Equivalence Test" begin
  line = (0, 20)
  ncells = 100

  mesh_uni = UniformGrid1D(line, ncells, :MEG6)
  mesh_curv = rectilinear_grid(line, ncells, :MEG6)

  # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.

  domain = mesh_uni.iterators.cell.domain

  # Cell center metrics
  @test isapprox(
    mesh_uni.cell_center_metrics.J[domain],
    mesh_curv.cell_center_metrics.J[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.x₁.ξ[domain],
    mesh_curv.cell_center_metrics.x₁.ξ[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.ξ.t[domain],
    mesh_curv.cell_center_metrics.ξ.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.ξ.x₁[domain],
    mesh_curv.cell_center_metrics.ξ.x₁[domain],
    atol=1e-10,
  )

  @test isapprox(
    mesh_uni.cell_center_metrics.ξ̂.t[domain],
    mesh_curv.cell_center_metrics.ξ̂.t[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.cell_center_metrics.ξ̂.x₁[domain],
    mesh_curv.cell_center_metrics.ξ̂.x₁[domain],
    atol=1e-10,
  )

  # Edge metrics
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ̂.x₁[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain],
    atol=1e-10,
  )
  @test isapprox(
    mesh_uni.edge_metrics.i₊½.ξ̂.t[domain],
    mesh_curv.edge_metrics.i₊½.ξ̂.t[domain],
    atol=1e-10,
  )
end
