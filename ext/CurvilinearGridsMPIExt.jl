module CurvilinearGridsMPIExt

using CurvilinearGrids
using CartesianDomains
using MPI
using KernelAbstractions

include("MPI/partitioned_rectlinear.jl")

export PartitionedRectlinearGrid

end