abstract type AbstractContinuousCurvilinearGrid1D <: AbstractCurvilinearGrid1D end
abstract type AbstractContinuousCurvilinearGrid3D <: AbstractCurvilinearGrid3D end
abstract type AbstractContinuousCurvilinearGrid2D <: AbstractCurvilinearGrid2D end

struct ContinuousCurvilinearGrid1D{T,A,B,C,EM,CM,E,BE,DBE,DS} <:
       AbstractContinuousCurvilinearGrid1D
  node_coordinates::A
  centroid_coordinates::A
  mapping_functions::B
  metric_functions_cache::C
  edge_metrics::EM
  cell_center_metrics::CM
  backend::BE
  diff_backend::DBE
  nhalo::Int
  discretization_scheme::DS
  discretization_scheme_name::Symbol
  iterators::E
end

struct ContinuousCurvilinearGrid2D{T,A,B,C,EM,CM,E,BE,DBE,DS} <:
       AbstractContinuousCurvilinearGrid2D
  node_coordinates::A
  centroid_coordinates::A
  mapping_functions::B
  metric_functions_cache::C
  edge_metrics::EM
  cell_center_metrics::CM
  backend::BE
  diff_backend::DBE
  nhalo::Int
  discretization_scheme::DS
  discretization_scheme_name::Symbol
  iterators::E
end

struct ContinuousCurvilinearGrid3D{T,A,B,C,EM,CM,E,BE,DBE,DS} <:
       AbstractContinuousCurvilinearGrid3D
  node_coordinates::A
  centroid_coordinates::A
  mapping_functions::B
  metric_functions_cache::C
  edge_metrics::EM
  cell_center_metrics::CM
  backend::BE
  diff_backend::DBE
  nhalo::Int
  discretization_scheme::DS
  discretization_scheme_name::Symbol
  iterators::E
end

include("1d.jl")
include("2d.jl")
include("3d.jl")
include("metric_cache.jl")
include("conserved_metrics.jl")
include("edge_interpolation.jl")
include("cell_center_derivs.jl")

function update_mapping_functions!(mesh, t, new_params, compute_metrics=true)
  mesh.mapping_function_params = new_params
  compute_node_coordinates!(mesh, t, new_params)
  compute_centroid_coordinates!(mesh, t, new_params)

  if compute_metrics
    compute_cell_metrics!(mesh, t, new_params)
    compute_edge_metrics!(mesh, t, new_params)
  end

  return nothing
end

function get_iterators(celldims::NTuple{N,Int}, nhalo::Int, global_cell_domain) where {N}
  cellCI = CartesianIndices(celldims .+ 2nhalo)
  nodeCI = CartesianIndices(celldims .+ 1 .+ 2nhalo)

  node = (full=nodeCI, domain=expand(nodeCI, -nhalo))
  cell = (full=cellCI, domain=expand(cellCI, -nhalo))

  if isnothing(global_cell_domain)
    global_domain = (node=node, cell=cell)
  else
    global_domain = (node=expand_upper(global_cell_domain, +1), cell=global_cell_domain)
  end
  return (; node, cell, nhalo, global_domain)
end