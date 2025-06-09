using Test, Random, CurvilinearGrids
using CurvilinearGrids.RectlinearArrays

@testset begin
    Random.seed!(1236)

    x = cumsum(rand(20));
    y = cumsum(rand(30));
    mesh_rect = rectlinear_grid(x, y, :MEG6)
    mesh_curv = rectlinear_grid(x, y, :MEG6, rectlinear=false)

    # Test to see if the domains of each array match. Note that the halo cell regions will NOT match, though this isn't an issue because the problem isn't defined there.
    
    domain = mesh_rect.iterators.cell.domain

    # Cell center metrics
    @test isapprox(mesh_rect.cell_center_metrics.J[domain], mesh_curv.cell_center_metrics.J[domain], atol=1e-14)

    @test isapprox(mesh_rect.cell_center_metrics.x₁.η[domain], mesh_curv.cell_center_metrics.x₁.η[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.x₁.ξ[domain], mesh_curv.cell_center_metrics.x₁.ξ[domain], atol=1e-14)

    @test isapprox(mesh_rect.cell_center_metrics.x₂.η[domain], mesh_curv.cell_center_metrics.x₂.η[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.x₂.ξ[domain], mesh_curv.cell_center_metrics.x₂.ξ[domain], atol=1e-14)

    @test isapprox(mesh_rect.cell_center_metrics.η.t[domain], mesh_curv.cell_center_metrics.η.t[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.η.x₁[domain], mesh_curv.cell_center_metrics.η.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.η.x₂[domain], mesh_curv.cell_center_metrics.η.x₂[domain], atol=1e-14)

    @test isapprox(mesh_rect.cell_center_metrics.η̂.t[domain], mesh_curv.cell_center_metrics.η̂.t[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.η̂.x₁[domain], mesh_curv.cell_center_metrics.η̂.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.η̂.x₂[domain], mesh_curv.cell_center_metrics.η̂.x₂[domain], atol=1e-14)

    @test isapprox(mesh_rect.cell_center_metrics.ξ.t[domain], mesh_curv.cell_center_metrics.ξ.t[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.ξ.x₁[domain], mesh_curv.cell_center_metrics.ξ.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.ξ.x₂[domain], mesh_curv.cell_center_metrics.ξ.x₂[domain], atol=1e-14)

    @test isapprox(mesh_rect.cell_center_metrics.ξ̂.t[domain], mesh_curv.cell_center_metrics.ξ̂.t[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.ξ̂.x₁[domain], mesh_curv.cell_center_metrics.ξ̂.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.cell_center_metrics.ξ̂.x₂[domain], mesh_curv.cell_center_metrics.ξ̂.x₂[domain], atol=1e-14)

    # Edge metrics
    @test isapprox(mesh_rect.edge_metrics.i₊½.ξ.x₁[domain], mesh_curv.edge_metrics.i₊½.ξ.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.ξ.x₂[domain], mesh_curv.edge_metrics.i₊½.ξ.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.η.x₁[domain], mesh_curv.edge_metrics.i₊½.η.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.η.x₂[domain], mesh_curv.edge_metrics.i₊½.η.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.ξ̂.x₁[domain], mesh_curv.edge_metrics.i₊½.ξ̂.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.ξ̂.x₂[domain], mesh_curv.edge_metrics.i₊½.ξ̂.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.ξ̂.t[domain], mesh_curv.edge_metrics.i₊½.ξ̂.t[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.η̂.x₁[domain], mesh_curv.edge_metrics.i₊½.η̂.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.η̂.x₂[domain], mesh_curv.edge_metrics.i₊½.η̂.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.i₊½.η̂.t[domain], mesh_curv.edge_metrics.i₊½.η̂.t[domain], atol=1e-14)

    @test isapprox(mesh_rect.edge_metrics.j₊½.ξ.x₁[domain], mesh_curv.edge_metrics.j₊½.ξ.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.ξ.x₂[domain], mesh_curv.edge_metrics.j₊½.ξ.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.η.x₁[domain], mesh_curv.edge_metrics.j₊½.η.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.η.x₂[domain], mesh_curv.edge_metrics.j₊½.η.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.ξ̂.x₁[domain], mesh_curv.edge_metrics.j₊½.ξ̂.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.ξ̂.x₂[domain], mesh_curv.edge_metrics.j₊½.ξ̂.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.ξ̂.t[domain], mesh_curv.edge_metrics.j₊½.ξ̂.t[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.η̂.x₁[domain], mesh_curv.edge_metrics.j₊½.η̂.x₁[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.η̂.x₂[domain], mesh_curv.edge_metrics.j₊½.η̂.x₂[domain], atol=1e-14)
    @test isapprox(mesh_rect.edge_metrics.j₊½.η̂.t[domain], mesh_curv.edge_metrics.j₊½.η̂.t[domain], atol=1e-14)
end
