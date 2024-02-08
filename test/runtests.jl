using CurvilinearGrids
using Test

@testset verbose = true "UnitTests" begin
  include("unit/test_indexing.jl")

  include("unit/test_meg6.jl")

  include("unit/test_2d.jl")

  include("unit/test_3d.jl")
end
