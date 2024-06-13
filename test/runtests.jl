using CurvilinearGrids
using Test
using StaticArrays
using BenchmarkTools
using LinearAlgebra
using StructArrays
using WriteVTK
using UnPack

@testset verbose = true "UnitTests" begin
  include("unit/test_indexing.jl")

  # include("unit/test_1d.jl")
  # include("unit/test_1d_axisymmetric.jl")

  include("unit/test_2d.jl")
  include("unit/test_2d_axisymmetric.jl")
  include("unit/perturb_mesh.jl")

  include("unit/test_3d.jl")
  include("unit/test_meg6.jl")
end
