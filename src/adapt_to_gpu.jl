using Adapt

function Adapt.adapt_structure(to, ds::MonotoneExplicit6thOrderDiscretization)
  xᵢ₊½ = Adapt.adapt_storage(to, ds.xᵢ₊½)
  ∂²x = Adapt.adapt_storage(to, ds.∂²x)
  ∂x = Adapt.adapt_storage(to, ds.∂x)
  metric = Adapt.adapt_storage(to, ds.metric)
  inner_deriv1 = Adapt.adapt_storage(to, ds.inner_deriv1)
  outer_deriv1 = Adapt.adapt_storage(to, ds.outer_deriv1)
  inner_deriv2 = Adapt.adapt_storage(to, ds.inner_deriv2)
  outer_deriv2 = Adapt.adapt_storage(to, ds.outer_deriv2)

  return MonotoneExplicit6thOrderDiscretization(
    xᵢ₊½, ∂²x, ∂x, metric, inner_deriv1, outer_deriv1, inner_deriv2, outer_deriv2
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
    grid.tiles,
    discretization_scheme,
    grid.onbc,
    grid.is_static,
    grid.is_orthogonal,
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
    grid.tiles,
    discretization_scheme,
    grid.onbc,
    grid.is_static,
    grid.is_orthogonal,
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
    grid.tiles,
    grid.snap_to_axis,
    grid.rotational_axis,
    grid.onbc,
    grid.is_static,
    grid.is_orthogonal,
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
    grid.tiles,
    discretization_scheme,
    grid.onbc,
    grid.is_static,
    grid.is_orthogonal,
  )
end
