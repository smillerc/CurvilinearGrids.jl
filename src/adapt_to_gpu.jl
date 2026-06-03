using Adapt
using KernelAbstractions

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

function Adapt.adapt_structure(to, grid::OrthogonalGrid)
  return _adapt_orthogonal_grid(
    KernelAbstractions.get_backend(grid.node_coordinates[1]), to, grid
  )
end

function _adapt_orthogonal_grid(::KernelAbstractions.GPU, to, grid::OrthogonalGrid)
  _ = to
  return grid
end

function _adapt_orthogonal_grid(::KernelAbstractions.Backend, to, grid::OrthogonalGrid)
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  cell_volumes = Adapt.adapt_structure(to, grid.cell_volumes)
  face_areas = Adapt.adapt_structure(to, grid.face_areas)

  return typeof(grid)(
    node_coordinates,
    centroid_coordinates,
    cell_volumes,
    grid.iterators,
    grid.domain_limits,
    face_areas,
    grid.nhalo,
  )
end

function Adapt.adapt_structure(to, cache::GridTypes.UnifiedMetricCache)
  data = Adapt.adapt_structure(to, cache.data)
  return GridTypes.UnifiedMetricCache(data, cache.valid, cache.mode)
end

function Adapt.adapt_structure(to, caches::GridTypes.UnifiedMetricCaches)
  cell = Adapt.adapt_structure(to, caches.cell)
  face = Adapt.adapt_structure(to, caches.face)
  return GridTypes.UnifiedMetricCaches(cell, face)
end

@inline _adapted_unified_backend(node_coordinates) =
  KernelAbstractions.get_backend(node_coordinates[1])

function Adapt.adapt_structure(to, grid::MappedGrid{N,T,CS,BT}) where {N,T,CS,BT}
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  face_coordinates = Adapt.adapt_structure(to, grid.face_coordinates)
  metric_caches = isnothing(grid.metric_caches) ? nothing : Adapt.adapt_structure(to, grid.metric_caches)
  backend = _adapted_unified_backend(node_coordinates)

  return MappedGrid{
    N,
    T,
    CS,
    BT,
    typeof(node_coordinates),
    typeof(centroid_coordinates),
    typeof(face_coordinates),
    typeof(grid.mapping_functions),
    typeof(grid.metric_functions_cache),
    typeof(backend),
    typeof(grid.diff_backend),
    typeof(grid.iterators),
    typeof(grid.state),
    typeof(metric_caches),
  }(
    node_coordinates,
    centroid_coordinates,
    face_coordinates,
    grid.mapping_functions,
    grid.metric_functions_cache,
    backend,
    grid.diff_backend,
    grid.nhalo,
    grid.iterators,
    grid.state,
    metric_caches,
  )
end

function Adapt.adapt_structure(to, grid::DiscreteGrid{N,T,CS,BT,IP}) where {N,T,CS,BT,IP}
  node_coordinates = Adapt.adapt_structure(to, grid.node_coordinates)
  centroid_coordinates = Adapt.adapt_structure(to, grid.centroid_coordinates)
  face_coordinates = Adapt.adapt_structure(to, grid.face_coordinates)
  metric_caches = isnothing(grid.metric_caches) ? nothing : Adapt.adapt_structure(to, grid.metric_caches)
  backend = _adapted_unified_backend(node_coordinates)

  return DiscreteGrid{
    N,
    T,
    CS,
    BT,
    IP,
    typeof(node_coordinates),
    typeof(centroid_coordinates),
    typeof(face_coordinates),
    typeof(grid.mapping_functions),
    typeof(grid.metric_functions_cache),
    typeof(backend),
    typeof(grid.diff_backend),
    typeof(grid.iterators),
    typeof(grid.state),
    typeof(metric_caches),
  }(
    node_coordinates,
    centroid_coordinates,
    face_coordinates,
    grid.mapping_functions,
    grid.metric_functions_cache,
    backend,
    grid.diff_backend,
    grid.nhalo,
    grid.iterators,
    grid.interpolation,
    grid.interpolants,
    grid.state,
    metric_caches,
  )
end
