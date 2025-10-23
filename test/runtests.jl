using BenchmarkTools
using CartesianDomains
using CurvilinearGrids
using CurvilinearGrids.RectilinearArrays
using DifferentiationInterface
using ForwardDiff: ForwardDiff
using KernelAbstractions
using LinearAlgebra
using StaticArrays
using StructArrays
using Test
using UnPack
using WriteVTK

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
  include("unit/test_continuous.jl")

  @info "Wall"
  include("unit/test_wall.jl")

  @info "Discretization Schemes"
  include("unit/test_discretizations.jl")
end
