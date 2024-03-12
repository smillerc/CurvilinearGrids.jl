using Adapt

function Adapt.adapt_structure(to, grid::CurvilinearGrid1D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  return CurvilinearGrid1D(
    node_coordinates,
    centroid_coordinates,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    grid.nhalo,
    grid.nnodes,
    grid.domain_limits,
    grid.iterators,
    grid._coordinate_funcs,
    grid._∂x∂ξ,
  )
end

function Adapt.adapt_structure(to, grid::CurvilinearGrid2D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  return CurvilinearGrid2D(
    node_coordinates,
    centroid_coordinates,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    grid.nhalo,
    grid.nnodes,
    grid.domain_limits,
    grid.iterators,
    grid._coordinate_funcs,
    grid._jacobian_matrix_func,
  )
end

function Adapt.adapt_structure(to, grid::CylindricalGrid2D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  edge_midpoint_coordinates = Adapt.adapt_structure(to, grid.edge_midpoint_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  return CylindricalGrid2D(
    node_coordinates,
    centroid_coordinates,
    edge_midpoint_coordinates,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    grid.nhalo,
    grid.nnodes,
    grid.domain_limits,
    grid.iterators,
    grid.snap_to_axis,
    grid._coordinate_funcs,
    grid._jacobian_matrix_func,
  )
end

function Adapt.adapt_structure(to, grid::CurvilinearGrid3D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  return CurvilinearGrid3D(
    node_coordinates,
    centroid_coordinates,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    grid.nhalo,
    grid.nnodes,
    grid.domain_limits,
    grid.iterators,
    grid.discretization_scheme,
    grid._coordinate_funcs,
    grid._jacobian_matrix_func,
  )
end
