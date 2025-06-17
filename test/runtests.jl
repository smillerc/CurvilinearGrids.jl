using CurvilinearGrids
using CurvilinearGrids.RectilinearArrays
using Test
using StaticArrays
using BenchmarkTools
using LinearAlgebra
using StructArrays
using WriteVTK
using UnPack

const test_gpu = false

include("unit/common.jl")
@testset verbose = true "UnitTests" begin
  include("unit/test_1d.jl")

  include("unit/test_2d.jl")
  include("unit/test_rectilinear_arrays.jl")
  include("unit/test_2d_axisymmetric.jl")
  include("unit/perturb_mesh.jl")

  # include("unit/test_3d.jl")
  include("unit/test_wall.jl")
  # include("unit/test_meg6.jl")
end
