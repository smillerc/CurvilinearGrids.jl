using CurvilinearGrids
using CurvilinearGrids.RectilinearArrays
using Test
using StaticArrays
using BenchmarkTools
using LinearAlgebra
using StructArrays
using WriteVTK
using UnPack
using KernelAbstractions

const test_gpu = false

include("unit/common.jl")

@testset verbose = true "UnitTests" begin
  @info "1D"
  include("unit/test_1d.jl")

  @info "2D"
  include("unit/test_2d.jl")

  @info "RectilinearArrays"
  include("unit/test_rectilinear_arrays.jl")

  @info "2D Axisymmetric"
  include("unit/test_2d_axisymmetric.jl")

  @info "Perturb example"
  include("unit/perturb_mesh.jl")

  @info "3D"
  include("unit/test_3d.jl")

  @info "Wall"
  include("unit/test_wall.jl")
  # include("unit/test_meg6.jl")

  @info "Read/Write .h5"
  include("unit/test_grid_read_write.jl")
end
