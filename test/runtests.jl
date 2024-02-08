using CurvilinearGrids
using Test
using StaticArrays
using BenchmarkTools

@testset verbose = true "UnitTests" begin
  include("unit/test_indexing.jl")

  include("unit/test_1d.jl")
  include("unit/test_2d.jl")
  include("unit/test_3d.jl")

  include("unit/test_meg6.jl")
end
