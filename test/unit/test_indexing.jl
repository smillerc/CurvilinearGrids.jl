@testset "Indexing functions" begin
  # include("../../src/indexing_utils.jl")
  using CurvilinearGrids.IndexingUtils
  domain = CartesianIndices((1:10, 4:8))

  @test lower_boundary_indices(domain, 1, +1) == CartesianIndices((2:2, 4:8))
  @test lower_boundary_indices(domain, 1, 0) == CartesianIndices((1:1, 4:8))
  @test lower_boundary_indices(domain, 1, -1) == CartesianIndices((0:0, 4:8))

  @test upper_boundary_indices(domain, 2, +1) == CartesianIndices((1:10, 9:9))
  @test upper_boundary_indices(domain, 2, 0) == CartesianIndices((1:10, 8:8))
  @test upper_boundary_indices(domain, 2, -1) == CartesianIndices((1:10, 7:7))

  @test expand(domain, 2, -1) == CartesianIndices((1:10, 5:7))
  @test expand(domain, 2, 0) == CartesianIndices((1:10, 4:8))
  @test expand(domain, 2, 2) == CartesianIndices((1:10, 2:10))

  @test expand(domain, 2) == CartesianIndices((-1:12, 2:10))

  @test expand_lower(domain, 2) == CartesianIndices((-1:10, 2:8))
  @test expand_upper(domain, 2) == CartesianIndices((1:12, 4:10))

  @test up(CartesianIndex((4, 5, 6)), 3, 2) == CartesianIndex(4, 5, 8)
  @test down(CartesianIndex((4, 5, 6)), 1, 2) == CartesianIndex(2, 5, 6)
  @test down(CartesianIndex((4, 5, 6)), 0, 2) == CartesianIndex(4, 5, 6)

  @test plus_minus(CartesianIndex((4, 5, 6)), 2, 2) == CartesianIndices((4:4, 3:7, 6:6))

  # non-uni form +/-: -1:+2
  @test plus_minus(CartesianIndex((4, 5, 6)), 2, (1, 2)) ==
    CartesianIndices((4:4, 4:7, 6:6))
  @test plus_minus(CartesianIndex((4, 5, 6)), 2, (1, 0)) ==
    CartesianIndices((4:4, 4:5, 6:6))

  # δ isn't exported...
  @test IndexingUtils.δ(1, CartesianIndex((4, 5, 6))) == CartesianIndex((1, 0, 0))
end
