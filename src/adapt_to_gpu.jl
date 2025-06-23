using Adapt

function Adapt.adapt_structure(to, grid::SphericalGrid1D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  discretization_scheme = Adapt.adapt_structure(to, grid.discretization_scheme)

  return SphericalGrid1D(
    node_coordinates,
    centroid_coordinates,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    grid.nhalo,
    grid.nnodes,
    grid.domain_limits,
    grid.iterators,
    discretization_scheme,
    grid.onbc,
    grid.snap_to_axis,
    grid.is_static,
    grid.discretization_scheme_name,
  )
end

function Adapt.adapt_structure(to, grid::CylindricalGrid1D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  discretization_scheme = Adapt.adapt_structure(to, grid.discretization_scheme)

  return CylindricalGrid1D(
    node_coordinates,
    centroid_coordinates,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    grid.nhalo,
    grid.nnodes,
    grid.domain_limits,
    grid.iterators,
    discretization_scheme,
    grid.onbc,
    grid.snap_to_axis,
    grid.is_static,
    grid.discretization_scheme_name,
  )
end

function Adapt.adapt_structure(to, grid::CurvilinearGrid1D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  discretization_scheme = Adapt.adapt_structure(to, grid.discretization_scheme)

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
    discretization_scheme,
    grid.is_static,
    grid.discretization_scheme_name,
  )
end

function Adapt.adapt_structure(to, grid::CurvilinearGrid2D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  discretization_scheme = Adapt.adapt_structure(to, grid.discretization_scheme)

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
    discretization_scheme,
    grid.is_static,
    grid.is_orthogonal,
    grid.discretization_scheme_name,
  )
end

function Adapt.adapt_structure(to, grid::AxisymmetricGrid2D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  edge_midpoint_coordinates = Adapt.adapt_structure(to, grid.edge_midpoint_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  return AxisymmetricGrid2D(
    node_coordinates,
    centroid_coordinates,
    edge_midpoint_coordinates,
    node_velocities,
    edge_metrics,
    cell_center_metrics,
    grid.nhalo,
    grid.nnodes,
    grid.iterators,
    grid.snap_to_axis,
    grid.rotational_axis,
    grid.onbc,
    grid.is_static,
    grid.is_orthogonal,
    grid.discretization_scheme_name,
  )
end

function Adapt.adapt_structure(to, grid::CurvilinearGrid3D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  discretization_scheme = Adapt.adapt_structure(to, grid.discretization_scheme)

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
    discretization_scheme,
    grid.is_static,
    grid.is_orthogonal,
    grid.discretization_scheme_name,
  )
end
