using CurvilinearGrids
using Makie
using Test

function _makie_rect_nodes_2d(x0, x1, y0, y1, ni, nj)
  xnodes = collect(range(x0, x1; length=ni + 1))
  ynodes = collect(range(y0, y1; length=nj + 1))
  x = [xnodes[i] for i in eachindex(xnodes), _ in eachindex(ynodes)]
  y = [ynodes[j] for _ in eachindex(xnodes), j in eachindex(ynodes)]
  return x, y
end

function _abutting_multiblock_2d(; nhalo=1)
  x1, y1 = _makie_rect_nodes_2d(0.0, 1.0, 0.0, 1.0, 4, 4)
  x2, y2 = _makie_rect_nodes_2d(1.0, 2.0, 0.0, 1.0, 4, 4)
  b1 = DiscreteGrid(x1, y1, nhalo)
  b2 = DiscreteGrid(x2, y2, nhalo)
  iface = (b1, :ihi) => (b2, :ilo)
  return MultiBlockMesh((b1, b2), (iface,); tolerance=1e-12)
end

@testset "Makie extension" begin
  ext = Base.get_extension(CurvilinearGrids, :CurvilinearGridsMakieExt)
  @test ext !== nothing

  @testset "Legacy 1D wireframe + mesh guard" begin
    r = collect(range(0.25, 2.0; length=24))
    grid = CylindricalGrid1D(r, :meg6, false)
    plt = Makie.wireframe(grid)
    @test plt isa Makie.FigureAxisPlot
    @test plt.axis isa Makie.Axis
    @test_throws ArgumentError Makie.mesh(grid)
  end

  @testset "Unified 2D wireframe and mesh" begin
    x, y = _makie_rect_nodes_2d(-1.0, 1.0, 0.0, 1.0, 8, 6)
    grid = DiscreteGrid(x, y, 1)
    wf = Makie.wireframe(grid)
    sf = Makie.mesh(grid)
    @test wf isa Makie.FigureAxisPlot
    @test sf isa Makie.FigureAxisPlot
    @test wf.axis isa Makie.Axis
    @test sf.axis isa Makie.Axis
  end

  @testset "Orthogonal 3D wireframe and mesh" begin
    r = collect(range(1.0, 1.5; length=6))
    θ = collect(range(0.3, 1.0; length=5))
    ϕ = collect(range(0.2, 1.1; length=5))
    grid = SphericalGrid3D(r, θ, ϕ, 0)
    wf = Makie.wireframe(grid)
    sf = Makie.mesh(grid)
    @test wf isa Makie.FigureAxisPlot
    @test sf isa Makie.FigureAxisPlot
    @test wf.axis isa Makie.Axis3
    @test sf.axis isa Makie.Axis3
  end

  @testset "MultiBlock color handling" begin
    mb = _abutting_multiblock_2d()
    wf = Makie.wireframe(mb; colors=(:red, :blue))
    sf = Makie.mesh(mb; colors=(:orange, :green))
    @test wf isa Makie.FigureAxisPlot
    @test sf isa Makie.FigureAxisPlot
    @test_throws ArgumentError Makie.wireframe(mb; colors=(:red,))
    @test_throws ArgumentError Makie.mesh(mb; colors=(:red,))
  end
end
