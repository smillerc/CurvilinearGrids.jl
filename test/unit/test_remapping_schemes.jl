using LinearAlgebra
using Test
using CurvilinearGrids

@testset "Legacy 1D remap helpers" begin
  source = CurvilinearGrid1D((0.0, 1.0), 8, :meg6)
  resized = change_resolution(source, 16)
  scaled = scale_resolution(source, 2.0)

  @test resized isa CurvilinearGrid1D
  @test scaled isa CurvilinearGrid1D
  @test length(coords(resized)) == 16
  @test length(coords(scaled)) == 17

  cell_data = collect(1.0:8.0)
  remapped = remap_cell_data(source, resized, cell_data)
  @test size(remapped) == size(resized.iterators.cell.domain)
  @test first(remapped) ≈ first(cell_data)
  @test last(remapped) ≈ last(cell_data)
end
