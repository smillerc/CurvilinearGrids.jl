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

include("unit/common.jl")

const test_gpu = false
const QUICK_TESTS = "--quick" in ARGS

@testset verbose = true "UnitTests" begin
  @info "1D"
  include("unit/test_1d.jl")
  include("unit/test_1d_continuous.jl")

  @info "2D"
  include("unit/test_2d.jl")
  include("unit/test_2d_continuous.jl")

  @info "RectilinearArrays"
  include("unit/test_rectilinear_arrays.jl")

  @info "2D Axisymmetric"
  include("unit/test_2d_axisymmetric.jl")

  @info "Orthogonal reduced grids"
  include("unit/test_orthogonal_grids.jl")

  @info "Perturb example"
  include("unit/perturb_mesh.jl")

  if QUICK_TESTS
    @info "Skipping long 3D grid tests (--quick)"
  else
    @info "3D"
    include("unit/test_3d.jl")
    include("unit/test_3d_continuous.jl")
  end

  @info "Wall"
  include("unit/test_wall.jl")

  @info "Discretization Schemes"
  include("unit/test_discretizations.jl")
  # include("unit/test_meg6.jl")

  @info "Unified Grid Types"
  include("unit/test_unified_grid_types.jl")
  include("unit/test_basis_transfer_api.jl")
  include("unit/test_unified_grid_backends.jl")
  include("unit/test_surface_grid.jl")

  @info "MultiBlock"
  include("unit/test_multiblock.jl")

  if !isnothing(Base.find_package("Makie"))
    @info "Makie Extension"
    include("unit/test_makie_ext.jl")
  else
    @info "Skipping Makie extension tests (Makie not available in test env)"
  end

  @info "Remapping"
  include("unit/test_remapping_schemes.jl")

  @info "Read/Write .h5"
  include("unit/test_grid_read_write.jl")
end
