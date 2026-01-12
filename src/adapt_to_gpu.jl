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
    grid.snap_to_axis,
    grid.is_static,
    grid.discretization_scheme_name,
    grid.halo_coords_included,
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
    grid.snap_to_axis,
    grid.is_static,
    grid.discretization_scheme_name,
    grid.halo_coords_included,
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
    grid.halo_coords_included,
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
    grid.halo_coords_included,
  )
end

function Adapt.adapt_structure(to, grid::AxisymmetricGrid2D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  node_velocities = Adapt.adapt_structure(to, grid.node_velocities)
  edge_metrics = Adapt.adapt_structure(to, grid.edge_metrics)
  cell_center_metrics = Adapt.adapt_structure(to, grid.cell_center_metrics)

  discretization_scheme = Adapt.adapt_structure(to, grid.discretization_scheme)

  return AxisymmetricGrid2D(
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
    grid.snap_to_axis,
    grid.rotational_axis,
    grid.is_static,
    grid.is_orthogonal,
    grid.discretization_scheme_name,
    grid.halo_coords_included,
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
    grid.halo_coords_included,
  )
end

function Adapt.adapt_structure(to, grid::SphericalGrid3D)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  cartesian_node_coordinates = Adapt.adapt_structure(to, grid.cartesian_node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  cell_volumes = Adapt.adapt_structure(to, grid.cell_volumes)
  face_areas = Adapt.adapt_structure(to, grid.face_areas)

  # cartesian_node_coords = (;
  #   x=KernelAbstractions.zeros(backend, T, nodedims),
  #   y=KernelAbstractions.zeros(backend, T, nodedims),
  #   z=KernelAbstractions.zeros(backend, T, nodedims),
  # )

  # spherical_centroid_coords = (;
  #   r=KernelAbstractions.zeros(backend, T, celldims[1]),
  #   θ=KernelAbstractions.zeros(backend, T, celldims[2]),
  #   ϕ=KernelAbstractions.zeros(backend, T, celldims[3]),
  # )

  # cell_volumes = KernelAbstractions.zeros(backend, T, celldims)

  # face_areas = (;
  #   i₊½=KernelAbstractions.zeros(backend, T, celldims),
  #   j₊½=KernelAbstractions.zeros(backend, T, celldims),
  #   k₊½=KernelAbstractions.zeros(backend, T, celldims),
  # )

  return SphericalGrid3D(
    node_coordinates,
    cartesian_node_coordinates,
    centroid_coordinates,
    cell_volumes,
    grid.iterators,
    grid.domain_limits,
    face_areas,
    grid.nhalo,
  )
end
